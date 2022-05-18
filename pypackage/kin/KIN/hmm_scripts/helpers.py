# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 11:33:59 2022

@author: divyaratan_popli
"""

import os
import numpy as np
import multiprocessing as mp

from .hmm_functions import *


def hmm_prep(targetfile, outfolder, datafolder, allrel, p_0):

    with open(targetfile) as f:
        libraries = [line.strip() for line in f]

    listf=[]
    for i, l1 in enumerate(libraries):
        for l2 in libraries[(i+1):]:
            s = '%s_._%s' %(str(l1).strip("['']"), str(l2).strip("['']"))
            listf.append(s)


    #make dirs

    if not(os.path.isdir(outfolder)):
        os.mkdir(outfolder)
    if not(os.path.isdir(outfolder+"likfiles/")):
        os.mkdir(outfolder+"likfiles/")
    if not(os.path.isdir(outfolder+"resfiles/")):
        os.mkdir(outfolder+"resfiles/")
    if not(os.path.isdir(outfolder+"gammafiles/")):
        os.mkdir(outfolder+"gammafiles/")


    if p_0==-9:
        pfile=datafolder + "hmm_parameters/p_0.txt"
        with open(pfile,"r") as f:
            p_0=float(f.read())

    dfile=datafolder + "input_diffs_hmm.csv"
    tfile=datafolder + "input_total_hmm.csv"

    p1thresh=datafolder+"input_hbd_hmm_total.csv"

    return libraries, listf, allrel, p_0, dfile, tfile, p1thresh



def hmm_run(hbdfolder, libraries, listf, instates, totalch, allrel, pfile, dfile, tfile, Afiles, p1thresh, outfolder, thresh, cores):

    pool = mp.Pool(cores)
    [pool.apply_async(hmm, args=(listind, hbdfolder, dfile, tfile,
                            np.array(listf), np.array(libraries), pfile,
                            Afiles, outfolder,
                            p1thresh, thresh, instates, allrel)) for listind in listf]
    pool.close()
    pool.join()

def hmm_results(outfolder, listf, allrel):

    likfile=[outfolder + "likfiles/{}.csv".format(n) for n in listf]
    relatable=getRelatable(filist=likfile, outfolder=outfolder, pairs=listf, rels=allrel)


    IBDadded=IBDstates(relatable=relatable, outfolder=outfolder)

    mergeRelatable(relf=IBDadded, outfolder=outfolder)



def hmm_all(targetfile, outfolder, datafolder, allrel, hbdfolder, thresh, cores, instates, totalch, p_0, Afiles):

    libraries, listf, allrel, p_0val, dfile, tfile, p1thresh = hmm_prep(targetfile = targetfile, outfolder = outfolder, datafolder = datafolder, allrel = allrel, p_0 = p_0)
    hmm_run(hbdfolder, libraries, listf, instates, totalch, allrel, p_0val, dfile, tfile, Afiles, p1thresh, outfolder, thresh, cores)
    hmm_results(outfolder, listf, allrel)
