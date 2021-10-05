# Test File Information

## Sources

## Extra Processing

```bash
# Cotter uvfits on 1254670392_avg with MWA flagging, both geometric and cable corrections
cotter \
  -m tests/data/1254670392_avg/1254670392.metafits \
  -o tests/data/1254670392_avg/1254670392.cotter.corrected.uvfits \
  -noantennapruning \
  -noflagautos \
  -noflagdcchannels \
  -nosbgains \
  -sbpassband tests/data/subband-passband-32ch-unitary.txt \
  -nostats \
  -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
  tests/data/1254670392_avg/1254670392*gpubox*.fits \
  | tee cotter-1254670392-uvfits-corrected.log
# Cotter measurement set on 1254670392_avg with MWA flagging, both geometric and cable corrections
cotter \
  -m tests/data/1254670392_avg/1254670392.metafits \
  -o tests/data/1254670392_avg/1254670392.cotter.corrected.ms \
  -noantennapruning \
  -noflagautos \
  -noflagdcchannels \
  -nosbgains \
  -sbpassband tests/data/subband-passband-32ch-unitary.txt \
  -nostats \
  -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
  tests/data/1254670392_avg/1254670392*gpubox*.fits \
  | tee cotter-1254670392-ms-corrected.log
# Cotter uvfits on 1254670392_avg with no flagging, or corrections
cotter \
  -m tests/data/1254670392_avg/1254670392.metafits \
  -o tests/data/1254670392_avg/1254670392.cotter.none.uvfits \
  -allowmissing \
  -norfi \
  -nostats \
  -nogeom \
  -noantennapruning \
  -nosbgains \
  -noflagautos \
  -noflagdcchannels \
  -nocablelength \
  -edgewidth 0 \
  -initflag 0 \
  -endflag 0 \
  -sbpassband tests/data/subband-passband-32ch-unitary.txt \
  -nostats \
  -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
  tests/data/1254670392_avg/1254670392*gpubox*.fits \
  | tee cotter-1254670392-uvfits-none.log
# Cotter measurement set on 1254670392_avg with no flagging, or corrections
cotter \
  -m tests/data/1254670392_avg/1254670392.metafits \
  -o tests/data/1254670392_avg/1254670392.cotter.none.ms \
  -allowmissing \
  -norfi \
  -nostats \
  -nogeom \
  -noantennapruning \
  -nosbgains \
  -noflagautos \
  -noflagdcchannels \
  -nocablelength \
  -edgewidth 0 \
  -initflag 0 \
  -endflag 0 \
  -sbpassband tests/data/subband-passband-32ch-unitary.txt \
  -nostats \
  -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
  tests/data/1254670392_avg/1254670392*gpubox*.fits \
  | tee cotter-1254670392-ms-none.log
```

then the following casa commands were used to create a truncated version of the table.

```python
tb.open('tests/data/1254670392_avg/1254670392.cotter.none.ms')
tb.copy('tests/data/1254670392_avg/1254670392.cotter.none.trunc.ms', norows=True)
tb.copyrows('tests/data/1254670392_avg/1254670392.cotter.none.trunc.ms', nrow=1)  
tb.close()
```