#!/usr/bin/env python

from os.path import abspath, exists, dirname
from os.path import join as path_join, exists as path_exists
from os import makedirs

from astropy.io import fits
from pprint import pformat
import numpy as np
import pandas as pd
import math
import re
import itertools
import sys
from argparse import ArgumentParser
import pyuvdata
from tabulate import tabulate

import numpy as np


def parse_args(argv):
    parser = ArgumentParser()
    parser.add_argument(
        "file",
        type=str
    )
    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    hdus = fits.open(args.file)
    print(f"-> hdus.info():")
    hdus.info()

    tile_header = hdus[1].data.dtype.names
    print(f"-> headers: {tile_header}")
    tile_data = hdus[1].data
    for row in tile_data:
        row[tile_header.index("Flag")] = 0
    hdus[1].data = tile_data
    print(tabulate(tile_data, headers=tile_header))
    hdus.writeto(args.file, overwrite=True)

if __name__ == '__main__':
    main(sys.argv[1:])
    # main(["tests/data/1254670392_avg/1254670392.metafits"])
