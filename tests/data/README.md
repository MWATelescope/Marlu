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
# Cotter measurement set on 1254670392_avg with no corrections, 4s, 80khz averaging
cotter \
  -m tests/data/1254670392_avg/1254670392.metafits \
  -o tests/data/1254670392_avg/1254670392.cotter.none.avg_4s_80khz.ms \
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
  -timeres 4 \
  -freqres 80 \
  tests/data/1254670392_avg/1254670392*gpubox*.fits \
  | tee cotter-1254670392-ms-none-avg_4s_80khz.log
```

Then, a version of the tables with only one baseline (truncated) is created with the following casa commands

```python
tb.open('tests/data/1254670392_avg/1254670392.cotter.none.ms')
tb.query('ANTENNA1 == 0 AND ANTENNA2 == 1 LIMIT 2').copy('tests/data/1254670392_avg/1254670392.cotter.none.trunc.ms', deep=True)
tb.close()
tb.open('tests/data/1254670392_avg/1254670392.cotter.none.avg_4s_80khz.ms')
tb.query('ANTENNA1 == 0 AND ANTENNA2 == 1 LIMIT 1').copy('tests/data/1254670392_avg/1254670392.cotter.none.avg_4s_80khz.trunc.ms', deep=True)
tb.close()
```

then the tables are dumped to csv with casa

```python
tb.open('tests/data/1254670392_avg/1254670392.cotter.none.trunc.ms')
exec(open('tests/data/casa_dump_ms.py').read())
```