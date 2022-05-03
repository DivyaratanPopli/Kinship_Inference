#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  1 22:15:29 2022

@author: divyaratan_popli
"""

import argparse
import multiprocessing as mp

from .KINgaroo_scripts import constants as C
from .KINgaroo_scripts import helpers as hel
from .KINgaroo_scripts import input_preparation_functions as inp

def cli():
    parser = argparse.ArgumentParser(description="Input generation pipeline for KIN")

    parser.add_argument('-bam', '--bamfiles_location',
                        type=str, metavar='',  required=True,
                        help='bamfiles directory')
    parser.add_argument('-bed', '--bedfile',
                        type=str, metavar='',  required=True,
                        help='path to bedfile')
    parser.add_argument('-i', '--interval',
                        type=str, metavar='', required=True,
                        help='Length of a genomic window in bases. Options:1e6, 1e7 (by default 1e7)')
    parser.add_argument('-T', '--target_location',
                        type=str, metavar='', required=True,
                        help='file with bamfile names that should be used without extension')
    parser.add_argument('-cnt', '--contam_parameter',
                        type=int, metavar='', required=True,
                        help='Enter 0 for no contamination correction. Enter 1 for contamination correction with divergence b/w target and contam. populations estimated from a vcf file (requires: contamination estimates file, indexed divergence vcf file, target individual, contaminating individual). Or enter divergence estimate s.t. 0<cnt<1')
    parser.add_argument('-c', '--cores',
                        type=int, metavar='',
                        help='Number of cores available')
    parser.add_argument('-t', '--threshold',
                        type=int, metavar='',
                        help='p_0 is estimated with all libraries that have atleast t number of informative windows (by default t=10)')
    
    parser.add_argument('-cest', '--contamination_estimates',
                        type=int, metavar='',
                        help='contamination estimates file with columns: name,contamination (%)')
    parser.add_argument('-vcf.gz', '--divergence_vcf.gz',
                        type=int, metavar='',
                        help='indexed compressed vcf file with an individual from target and contaminating populations each. Diploid Genotypes (GT) should be represented')
    parser.add_argument('-tar', '--target_ind',
                        type=int, metavar='',
                        help='Name of target individual in divergence_vcf file')
    parser.add_argument('-cont', '--contaminating_ind',
                        type=int, metavar='',
                        help='Name of contaminating individual in divergence_vcf file')
    parser.add_argument('-r', '--roh',
                        type=int, metavar='',
                        help='Enter 1 if you need ROH estimates. Enter 0 if you already have the positions of ROH tracts (by default:1).')
    return parser.parse_args()


def main():
    args = cli()

    if args.cores is None:
        cores = mp.cpu_count()
    else:
        cores = args.cores
    if args.threshold is None:
        thresh = 10
    else:
        thresh = args.threshold
    if args.interval is None:
        interval = 1e7
    else:
        interval = args.interval
    if args.roh is None:
        roh = 1
    else:
        roh = args.roh
    libraries, listf, dwins, twins, id_dwins, id_twins, chrmlist = hel.pipeline1(targetsfile = args.target_location, 
                                                                                 bedfile = args.bedfile, 
                                                                                 cores = args.cores, 
                                                                                 rawbams = args.bamfiles_location, 
                                                                                 interval = interval)
    
    if args.contam_parameter==0:
            diff_cor, total_cor, id_diff_cor, id_total_cor = dwins, twins, id_dwins, id_twins
    elif args.contam_parameter!=0:
        inp.contamFile(infile=args.contamination_estimates, outfile=C.contam_est, targets=libraries, idfile=C.id_contam_est)
        if args.contam_parameter==1:
            div= hel.divergence_vcf(args.divergence_vcf.gz, args.bedfile, args.target_ind, args.contaminating_ind)
        elif args.contam_parameter>0 and args.contam_parameter<1:
            div = args.contam_parameter
        
    diff_cor,total_cor = inp.contamAll(diff=dwins, total=twins, cfile=C.contam_est, p_c=div)
    id_diff_cor, id_total_cor = inp.contamAll(diff=id_dwins, total=id_twins, cfile=C.id_contam_est, p_c=div)

    diff_cor, total_cor, id_diff_cor, id_total_cor, p1 = hel.data2p(divergence=div, diff_cor=diff_cor, 
                                                            total_cor=total_cor, id_diff_cor=id_diff_cor, 
                                                            id_total_cor=id_total_cor, libraries=libraries, listf=listf, 
                                                            hmm_param=C.hmm_param, thresh=thresh, outdiff=C.outdiff, outtotal=C.outtotal, 
                                                            id_outdiff=C.outdiff, id_outtotal=C.outtotal)

    if roh==1:      
        hel.run_hmm(diff=diff_cor, total=total_cor, chrm1=chrmlist, p1=p1, hbdf=C.hbdf, likf=C.likf, libraries=libraries, cores=cores)
    


if __name__ == '__main__':
    main()
