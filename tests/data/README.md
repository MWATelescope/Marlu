# Test File Information

## Sources

## Clear antenna flags

the obs we're using has all antennas flagged, this unflags them.

```bash
python3 tests/data/clear_ant_flags.py tests/data/1254670392_avg/1254670392.metafits
```

## Run through Cotter

```bash
# Cotter measurement set on 1254670392_avg with no corrections
cotter \
  -m tests/data/1254670392_avg/1254670392.metafits \
  -o tests/data/1254670392_avg/1254670392.cotter.none.ms \
  -allowmissing \
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
  -flag-strategy /usr/share/aoflagger/strategies/mwa-default.lua \
  tests/data/1254670392_avg/1254670392*gpubox*.fits \
  | tee cotter-1254670392-ms-none.log
```

then the following casa commands were used to create a truncated version of the table.

```python
tb.open('tests/data/1254670392_avg/1254670392.cotter.none.ms')
tb.copy('tests/data/1254670392_avg/1254670392.cotter.none.trunc.ms', norows=True)
tb.copyrows('tests/data/1254670392_avg/1254670392.cotter.none.trunc.ms', startrowin=1, nrow=1)  
tb.close()
```

