# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 11:33:59 2022

@author: divyaratan_popli
"""

import os
import numpy as np
import multiprocessing as mp

import constants as C
import hmm_functions as hf


def hmm_prep(targetfile, outfolder, datafolder, paramfolder, allrel):
    
    with open(targetfile) as f:
        libraries = [line.strip() for line in f]
    
    listf=[]
    for i, l1 in enumerate(libraries):
        for l2 in libraries[(i+1):]:
            s = '%s_._%s' %(str(l1).strip("['']"), str(l2).strip("['']"))
            listf.append(s)

    instates=C.STATES
    totalch=C.CHRM
    
    allrel=C.RELS
    #make dirs
    
    if not(os.path.isdir(outfolder)):
        os.mkdir(outfolder)
    if not(os.path.isdir(outfolder+"likfiles/")):
        os.mkdir(outfolder+"likfiles/")
    if not(os.path.isdir(outfolder+"resfiles/")):
        os.mkdir(outfolder+"resfiles/")
    if not(os.path.isdir(outfolder+"gammafiles/")):
        os.mkdir(outfolder+"gammafiles/")



    pfile=datafolder + "hmm_parameters_fil0/p1_file"
    dfile=datafolder + "mergedwin_remHighdiv_fil0/pw_diff.csv"
    tfile=datafolder + "mergedwin_remHighdiv_fil0/pw_total.csv"
    
    Afiles=[paramfolder + "transition_matrix_{}.csv".format(n) for n in allrel]
    p1thresh=datafolder+"identicalmergedwin_remHighdiv_fil0/id_total.csv"
    
    return libraries, listf, instates, totalch, allrel, pfile, dfile, tfile, Afiles, p1thresh



def hmm_run(hbdfolder, libraries, listf, instates, totalch, allrel, pfile, dfile, tfile, Afiles, p1thresh, outfolder, thresh, cores):
    pool = mp.Pool(cores)
    [pool.apply_async(hf.hmm, args=(listind, hbdfolder, dfile, tfile,
                            np.array(listf), np.array(libraries), pfile,
                            Afiles, outfolder,
                            p1thresh, thresh, instates, allrel)) for listind in listf]
    pool.close()
    pool.join()

def hmm_results(outfolder, listf, allrel):

    likfile=[outfolder + "likfiles/{}".format(n) for n in listf]
    relatable=hf.getRelatable(filist=likfile, outfolder=outfolder, pairs=listf, rels=allrel)


    IBDadded=hf.IBDstates(relatable=relatable, outfolder=outfolder)

    hf.mergeRelatable(relf=IBDadded, outfolder=outfolder)



def hmm_all(targetfile, outfolder, datafolder, paramfolder, allrel, hbdfolder, thresh, cores):

    libraries, listf, instates, totalch, allrel, pfile, dfile, tfile, Afiles, p1thresh = hmm_prep(targetfile = targetfile, outfolder = outfolder, datafolder = datafolder, paramfolder = paramfolder, allrel = allrel)
    hmm_run(hbdfolder, libraries, listf, instates, totalch, allrel, pfile, dfile, tfile, Afiles, p1thresh, outfolder, thresh, cores)
    hmm_results(outfolder, listf, allrel)
