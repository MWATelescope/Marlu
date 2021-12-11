//! Create a measurement set, and write synthetic visiblity data into the DATA column.

use marlu::{io::ms::MeasurementSetWriter, RADec, Complex};
use rubbl_casatables::{Array, Table, TableOpenMode};
use itertools::iproduct;

use lazy_static::lazy_static;

lazy_static! {
    static ref DEFAULT_TABLES_GZ: &'static [u8] = include_bytes!("../data/default_tables.tar.gz");
}

fn main() {
    // DELETEME:
    coz::scope!("main");

    // Set the shape of the synthetic data
    let num_ants = 128;
    let num_baselines = num_ants * (num_ants - 1) / 2;
    let num_timesteps = 10;
    let num_rows = num_baselines * num_timesteps;
    let num_channels = 24 * 32;
    let num_pols = 4;
    let data_shape = (num_channels, num_pols);

    // Create a new temporary directory to write to each time
    let tmp_dir = tempfile::tempdir().unwrap();
    let table_path = tmp_dir.path().join("table.ms");

    // set up an ms writer
    let phase_centre = RADec::new(0., 0.);
    let ms_writer = MeasurementSetWriter::new(table_path.clone(), phase_centre, None);
    ms_writer.decompress_default_tables().unwrap();
    ms_writer.decompress_source_table().unwrap();
    ms_writer.add_cotter_mods(num_channels);

    // Create an array to store synthetic visiblity data, re-used each row.
    let mut uvw_tmp = vec![0.; 3];
    let sigma_tmp = vec![1.; 4];
    let mut data_tmp = Array::<Complex<f32>, _>::zeros(data_shape);
    let mut weights_tmp = Array::<f32, _>::zeros(data_shape);
    let mut flags_tmp = Array::from_elem(data_shape, false);

    let mut table = Table::open(table_path.clone(), TableOpenMode::ReadWrite).unwrap();

    {
        coz::scope!("add rows");
        table.add_rows(num_rows).unwrap();
    }

    // Write synthetic visibility data to the data column of the table.
    for (row_index, (timestep_index, baseline_index)) in iproduct!(0..num_timesteps, 0..num_baselines).enumerate() {
        // DELETEME:
        coz::scope!("row loop");

        // Calculate the uvw coordinates for this row.
        uvw_tmp[0] = row_index as _;
        uvw_tmp[1] = baseline_index as _;
        uvw_tmp[2] = timestep_index as _;

        // Calculate the weights for this row.
        weights_tmp.column_mut(0).fill(row_index as _);
        weights_tmp.column_mut(1).fill(baseline_index as _);
        weights_tmp.column_mut(2).fill(timestep_index as _);
        weights_tmp.column_mut(3).fill(1.);

        // each element in the visibility array is a complex number, whose real component is
        // the row index and whose imaginary component is the element index.
        data_tmp.iter_mut().enumerate().for_each(|(i, v)| {
            *v = Complex::new(row_index as f32, i as f32);
        });
        flags_tmp.iter_mut().enumerate().for_each(|(i, v)| {
            *v = i % 2 == 0 as _;
        });

        let antenna1 = baseline_index % num_ants;
        let antenna2 = baseline_index / num_ants;

        // Write the row to the table.
        ms_writer.write_main_row(
            &mut table,
            row_index as _,
            timestep_index as _,
            timestep_index as _,
            antenna1 as _,
            antenna2 as _,
            0, 
            &uvw_tmp,
            1.,
            -1,
            1,
            -1,
            &sigma_tmp,
            &data_tmp,
            &flags_tmp,
            &weights_tmp,
            row_index % 2 == 0 as _,
        ).unwrap();
    }
}
