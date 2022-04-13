// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Convenience code for interfacing with CUDA.

// This module should only be accessed by lib.rs if the "cuda" feature is
// enabled, so we don't need to have conditional compilation here.

use std::ffi::{c_void, CStr, CString};

use thiserror::Error;

/// The length of the error strings used to get error messages from CUDA
/// functions.
pub const ERROR_STR_LENGTH: usize = 1024;

/// A Rust-managed pointer to CUDA device memory. When this is dropped,
/// [`cuda_runtime_sys::cudaFree`] is called on the pointer.
#[derive(Debug)]
pub struct DevicePointer<T> {
    ptr: *mut T,

    /// The number of elements of `T` allocated against the pointer.
    ///
    /// Yeah, I know a pointer isn't an array.
    num_elements: usize,
}

impl<T> DevicePointer<T> {
    /// Allocate a number of bytes on the device.
    ///
    /// # Safety
    ///
    /// This function interfaces directly with the CUDA API. Rust errors attempt
    /// to catch problems but there are no guarantees.
    pub unsafe fn malloc(size: usize) -> Result<DevicePointer<T>, CudaError> {
        let mut d_ptr = std::ptr::null_mut();
        cuda_runtime_sys::cudaMalloc(&mut d_ptr, size);
        peek_and_sync(CudaCall::Malloc)?;
        Ok(Self {
            ptr: d_ptr.cast(),
            num_elements: size / std::mem::size_of::<T>(),
        })
    }

    /// Get the number of elements of `T` that have been allocated on the
    /// device. The number of bytes allocated is `num_elements *
    /// std::mem::size_of::<T>()`.
    pub fn get_num_elements(&self) -> usize {
        self.num_elements
    }

    /// Copy a slice of data to the device. Any type is allowed, and the returned
    /// pointer is to the device memory.
    ///
    /// # Safety
    ///
    /// This function interfaces directly with the CUDA API. Rust errors attempt
    /// to catch problems but there are no guarantees.
    pub unsafe fn copy_to_device(v: &[T]) -> Result<DevicePointer<T>, CudaError> {
        let size = v.len() * std::mem::size_of::<T>();
        let d_ptr = Self::malloc(size)?;
        cuda_runtime_sys::cudaMemcpy(
            d_ptr.get_mut() as *mut c_void,
            v.as_ptr().cast(),
            size,
            cuda_runtime_sys::cudaMemcpyKind::cudaMemcpyHostToDevice,
        );
        peek_and_sync(CudaCall::CopyToDevice)?;
        Ok(d_ptr)
    }

    /// Copy a slice of data from the device. The amount of data copied depends
    /// on the length of `v`, so if in doubt, the size of the device allocation
    /// should be checked with `DevicePointer::get_num_elements`.
    ///
    /// # Safety
    ///
    /// This function interfaces directly with the CUDA API. Rust errors attempt
    /// to catch problems but there are no guarantees.
    pub unsafe fn copy_from_device(&self, v: &mut [T]) -> Result<(), CudaError> {
        let size = v.len() * std::mem::size_of::<T>();
        cuda_runtime_sys::cudaMemcpy(
            v.as_mut_ptr().cast(),
            self.ptr.cast(),
            size,
            cuda_runtime_sys::cudaMemcpyKind::cudaMemcpyDeviceToHost,
        );
        peek_and_sync(CudaCall::CopyFromDevice)
    }

    /// Overwrite the device memory allocated against this [`DevicePointer`]
    /// with new memory. The amount of memory `v` must match what is allocated
    /// on against this [`DevicePointer`].
    ///
    /// # Safety
    ///
    /// This function interfaces directly with the CUDA API. Rust errors attempt
    /// to catch problems but there are no guarantees.
    pub unsafe fn overwrite(&mut self, v: &[T]) -> Result<(), CudaError> {
        if v.len() != self.num_elements {
            return Err(CudaError::SizeMismatch);
        }
        let size = v.len() * std::mem::size_of::<T>();
        cuda_runtime_sys::cudaMemcpy(
            self.get_mut() as *mut c_void,
            v.as_ptr().cast(),
            size,
            cuda_runtime_sys::cudaMemcpyKind::cudaMemcpyHostToDevice,
        );
        peek_and_sync(CudaCall::CopyToDevice)
    }

    /// Get a const pointer to the device memory.
    pub fn get(&self) -> *const T {
        self.ptr as *const T
    }

    /// Get a mutable pointer to the device memory.
    pub fn get_mut(&self) -> *mut T {
        self.ptr
    }
}

impl<T> Drop for DevicePointer<T> {
    fn drop(&mut self) {
        unsafe {
            cuda_runtime_sys::cudaFree(self.ptr.cast());
        }
    }
}

#[derive(Error, Debug)]
pub enum CudaError {
    #[error("When overwriting, the new amount of memory did not equal the old amount")]
    SizeMismatch,

    #[error("cudaMemcpy to device failed: {0}")]
    CopyToDevice(String),

    #[error("cudaMemcpy from device failed: {0}")]
    CopyFromDevice(String),

    #[error("cudaMalloc error: {0}")]
    Malloc(String),

    #[error("CUDA kernel error: {0}")]
    Kernel(String),
}

/// Turn a non-zero exit code and an error string pointer (which was originally
/// allocated by Rust as a [`CString`]) into a Result. The error can only be the
/// `CudaError::Kernel` variant. The memory associated with the `error_str`
/// pointer is consumed by this function.
///
/// # Safety
///
/// This function assumes that `error_str` is a valid C string. If `error_str`
/// is used again after calling this function, undefined behaviour is expected.
pub unsafe fn cuda_status_to_error(
    exit_code: i32,
    error_str: *mut std::os::raw::c_char,
) -> Result<(), CudaError> {
    // This stuff to convert the pointer into a CString is necessary to properly
    // deallocate the Rust-allocated string.
    let error_str = CString::from_raw(error_str)
        .into_string()
        // This assumes that the foreign code doesn't corrupt the string.
        .unwrap();
    if exit_code == 0 {
        Ok(())
    } else {
        Err(CudaError::Kernel(error_str))
    }
}

#[derive(Clone, Copy)]
pub enum CudaCall {
    Malloc,
    CopyToDevice,
    CopyFromDevice,
}

/// Run [`cuda_runtime_sys::cudaPeekAtLastError`] and
/// [`cuda_runtime_sys::cudaDeviceSynchronize`]. If either of these calls return
/// an error, it is converted to a Rust error and returned from this function.
/// The single argument describes what the just-performed operation was and
/// makes the returned error a helpful one.
///
/// # Safety
///
/// This function interfaces directly with the CUDA API. Rust errors attempt to
/// catch problems but there are no guarantees.
pub unsafe fn peek_and_sync(cuda_call: CudaCall) -> Result<(), CudaError> {
    let code = cuda_runtime_sys::cudaPeekAtLastError();
    if code != cuda_runtime_sys::cudaError::cudaSuccess {
        let c_str = CStr::from_ptr(cuda_runtime_sys::cudaGetErrorString(code));
        let s = c_str.to_str().unwrap().to_string();
        return Err(match cuda_call {
            CudaCall::Malloc => CudaError::Malloc(s),
            CudaCall::CopyToDevice => CudaError::CopyToDevice(s),
            CudaCall::CopyFromDevice => CudaError::CopyFromDevice(s),
        });
    }

    let code = cuda_runtime_sys::cudaDeviceSynchronize();
    if code != cuda_runtime_sys::cudaError::cudaSuccess {
        let c_str = CStr::from_ptr(cuda_runtime_sys::cudaGetErrorString(code));
        let s = c_str.to_str().unwrap().to_string();
        return Err(match cuda_call {
            CudaCall::Malloc => CudaError::Malloc(s),
            CudaCall::CopyToDevice => CudaError::CopyToDevice(s),
            CudaCall::CopyFromDevice => CudaError::CopyFromDevice(s),
        });
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use ndarray::Array1;
    use serial_test::serial;

    #[test]
    fn copy_to_and_from_device_succeeds() {
        unsafe {
            let heap = vec![0_u32; 100];
            let result = DevicePointer::copy_to_device(&heap);
            assert!(result.is_ok(), "Couldn't copy data to device memory");
            let d_ptr = result.unwrap();
            let mut heap2 = vec![1_u32; d_ptr.get_num_elements()];
            let result = d_ptr.copy_from_device(&mut heap2);
            assert!(result.is_ok(), "Couldn't copy data from device memory");
            result.unwrap();
            assert_abs_diff_eq!(Array1::from(heap), Array1::from(heap2));

            const LEN: usize = 100;
            let stack = [0_u32; LEN];
            let result = DevicePointer::copy_to_device(&stack);
            assert!(result.is_ok(), "Couldn't copy data to device memory");
            let d_ptr = result.unwrap();
            let mut stack2 = [1_u32; LEN];
            let result = d_ptr.copy_from_device(&mut stack2);
            assert!(result.is_ok(), "Couldn't copy data from device memory");
            result.unwrap();
            assert_abs_diff_eq!(Array1::from(stack.to_vec()), Array1::from(stack2.to_vec()));

            // Verify the copy_from_behaviour behaviour that the number of
            // elements copied depends on the host's array length.
            let mut stack3 = [1_u32; 100];
            let result = d_ptr.copy_from_device(&mut stack3[..50]);
            assert!(result.is_ok(), "Couldn't copy data from device memory");
            result.unwrap();
            assert_abs_diff_eq!(Array1::from(stack3[..50].to_vec()), Array1::zeros(50));
            assert_abs_diff_eq!(Array1::from(stack3[50..].to_vec()), Array1::ones(50));
        }
    }

    #[test]
    #[serial]
    fn cuda_malloc_huge_fails() {
        let size = 1024_usize.pow(4); // 1 TB;
        let result: Result<DevicePointer<u8>, CudaError> = unsafe { DevicePointer::malloc(size) };
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert_eq!(err.to_string(), "cudaMalloc error: out of memory");
    }

    #[test]
    fn copy_from_non_existent_pointer_fails() {
        let d_ptr: DevicePointer<u8> = DevicePointer {
            ptr: std::ptr::null_mut::<u8>(),
            num_elements: 1,
        };
        let mut dest = [0; 100];
        let result = unsafe { d_ptr.copy_from_device(&mut dest) };
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert_eq!(
            err.to_string(),
            "cudaMemcpy from device failed: invalid argument"
        );
    }
}
