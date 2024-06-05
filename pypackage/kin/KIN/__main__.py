#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 17:52:45 2022

@author: divyaratan_popli
"""

import argparse
import multiprocessing as mp

from .hmm_scripts import constants as C
from .hmm_scripts import helpers


def cli():
    parser = argparse.ArgumentParser(description="Relatedness and IBD estimates")

    parser.add_argument('-I', '--input_location',
                        type=str, required=True,
                        help='input files location')
    parser.add_argument('-O', '--output_location',
                        type=str, required=True,
                        help='Output files location')
    parser.add_argument('-r', '--ROH_file_location',
                        type=str,
                        help='ROH files location')
    parser.add_argument('-c', '--cores',
                        type=int,
                        help='Number of cores available')
    parser.add_argument('-t', '--threshold',
                        type=int,
                        help='Minimum number of sites in a window for ROH implementation')
    parser.add_argument('-p', '--diversity_parameter_p_0',
                        type=float, metavar='p_0',
                        help='Input p_0 parameter, if you do not want to calculate it from given samples (Keep it same as that for KINgaroo)')

    return parser.parse_args()


def main():
    args = cli()

    if args.ROH_file_location is None:
        hbdfolder = args.input_location
    else:
        hbdfolder = args.ROH_file_location

    if args.threshold is None:
        thresh = 10
    else:
        thresh = args.threshold

    if args.cores is None:
        cores = mp.cpu_count()
    else:
        cores = args.cores

    if args.diversity_parameter_p_0 is None:
        p_0 = -9
    else:
        p_0 = args.diversity_parameter_p_0


    with open(args.input_location+"/interval.txt","r") as f:
        interval1=int(f.read())


    if interval1==10000000:
        A_interval=C.A_ALL10
    elif interval1==1000000:
        A_interval=C.A_ALL1

    helpers.hmm_all(
        targetfile=args.input_location + "/target_samples.txt",
        outfolder=args.output_location + "/",
        datafolder=args.input_location + "/",
        allrel=C.RELS,
        hbdfolder=hbdfolder+ "/",
        thresh=thresh,
        cores=cores,
        instates=C.STATES,
        p_0=p_0,
        Afiles=A_interval
    )


if __name__ == '__main__':
    main()
