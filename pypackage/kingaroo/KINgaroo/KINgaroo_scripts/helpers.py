#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 11:49:28 2022

@author: divyaratan_popli
"""
import pysam
import pandas as pd
import numpy as np
import os
import subprocess
import multiprocessing as mp
from .input_preparation_functions import *
from .hbd_hmm_functions import *

def test_input(bedfile, rawbams, targetsfile):
    print("testing the input files...")

    bed=pd.read_csv(bedfile, sep="\t", header=None, index_col=False, low_memory=False).loc[:10,:]
    with pd.option_context('display.max_rows', len(bed.index), 'display.max_columns', len(bed.columns)):
        bed.to_csv("test_bed10.bed", sep='\t', header=None, index=False)

    with open(targetsfile) as f:
        lib = [line.strip() for line in f]

    np.savetxt(fname="test_targets.txt", X=lib, delimiter=',', fmt="%s")

    bamf=rawbams + lib[0] + '.bam'
    command = 'samtools view -b %s %s -o %s' %(bamf, 22, "test_bam22.bam")
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    print("Input files are OK.")
    print("Please further check that these files are not empty: test_targets.txt (names of all libraries), test_bam22.bam (bamfile for first library with chromosome 22), test_bed10.bed (first 10 rows of bedfile)")
#creating initial indexes:


def prep_function(targetsfile, splitbams, bedfiles, hapProbs, hmm_param, hbd, lik):

    with open(targetsfile) as f:
        libraries = [line.strip() for line in f]

    listf=[]
    for i, l1 in enumerate(libraries):
        for l2 in libraries[(i+1):]:
            s = '%s_%s' %(str(l1).strip("['']"), str(l2).strip("['']"))
            listf.append(s)


    if not(os.path.isdir(splitbams)):
        os.mkdir(splitbams)
    if not(os.path.isdir(bedfiles)):
        os.mkdir(bedfiles)
    if not(os.path.isdir(hapProbs)):
        os.mkdir(hapProbs)
    if not(os.path.isdir(hmm_param)):
        os.mkdir(hmm_param)
    if not(os.path.isdir(hbd)):
        os.mkdir(hbd)
    if not(os.path.isdir(lik)):
        os.mkdir(lik)

    np.savetxt(fname="target_samples.txt", X=libraries, delimiter="\t", fmt='%s')

    return libraries, listf

def split_bed(bedfile, totalch, bedfiles):
    bed=pd.read_csv(bedfile, sep="\t", header=None, index_col=False, low_memory=False)
    for b in totalch:
        chrm_bed=bed.loc[bed[0].astype(str)==str(b),:]

        bed_chrmf=bedfiles+'bedfile_chrm%s.bed' %(str(b))
        with pd.option_context('display.max_rows', len(chrm_bed.index), 'display.max_columns', len(chrm_bed.columns)):
            chrm_bed.to_csv(bed_chrmf, sep='\t', header=None, index=False)


#creating initial indexes:
def init_index(rawbams, lib):
    bamf=rawbams + lib + '.bam'
    print("Indexing %s...\n" %(lib))
    pysam.index(bamf)

def create_lib_chrm(libraries,totalch):
    lib_chrm_all=[]
    for lib in libraries:
        for ch in totalch:
            lib_chrm_all.append((lib,ch))
    return lib_chrm_all


def bam_filter(rawbams, splitbams, bedfiles, hapProbs, lib_chrm, sort1):

    lib=lib_chrm[0]
    chrm=int(lib_chrm[1])

    if str(sort1)=='1':
        print("sorting library %s chromosome %s...\n" %(lib,chrm))
        bamsplit_chr=splitbams+'%s_chrm%s.bam' %(lib,chrm)
    elif str(sort1)=='0':
        bamsplit_chr=splitbams+'%s_chrm%s.sorted.bam' %(lib,chrm)
    bamf=rawbams + lib + '.bam'
    command = 'samtools view -b %s %s -o %s' %(bamf, chrm, bamsplit_chr)
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

    bamsorted=splitbams+'%s_chrm%s.sorted.bam' %(lib,chrm)
    if str(sort1)=='1':
        pysam.sort("-o", bamsorted, bamsplit_chr)
        pysam.index(bamsorted)
    elif str(sort1)=='0':
        pysam.index(bamsplit_chr)

    #input file generation
    bedfile=bedfiles + 'bedfile_chrm%s.bed' %(chrm)
    baifile=bamsorted+'.bai'

    outprob=hapProbs + 'hapProbs_%s_chrm%s_probs.csv' %(lib,chrm)
    outdiff=hapProbs + 'hapProbs_%s_chrm%s_diffs.csv' %(lib,chrm)

    makeHapProbs(inbam=bamsorted, inbai=baifile, inbed=bedfile, outprob=outprob, outdiff_id=outdiff)

#getting splitbams for all libs and chrm in parallel
def parallel_indexes(rawbams, libraries, cores):
    pool = mp.Pool(cores)
    [pool.apply_async(init_index, args=(rawbams, lib)) for lib in libraries]
    pool.close()
    pool.join()


def parallel_bamfilter(rawbams, splitbams, bedfiles, hapProbs, lib_chrm_all, cores, sort1):
    pool = mp.Pool(cores)
    [pool.apply_async(bam_filter, args=(rawbams, splitbams, bedfiles, hapProbs, lib_chrm, sort1)) for lib_chrm in lib_chrm_all]
    pool.close()
    pool.join()
    print("Finished indexing, sorting and splitting bamfiles...")

def get_merged_chrm(libraries, chrm, interval):
    print("Creating input files from chromosome %s..." %(chrm))

    error_in="Can not process input data. Please check all input files. In particular, bedfile and bamfiles should have chromosomes listed as 1,2,... and the bedfile should be tab-separated. Please check the README.md for example input files."
    try:
        pslist=[("hapProbs/hapProbs_{}_chrm%s_probs.csv" %(chrm)).format(n) for n in libraries]
        probs_list, pos_list, chrm_list = hapProbsAll(haplist=pslist, hap='noid')

        id_diffs =[("hapProbs/hapProbs_{}_chrm%s_diffs.csv" %(chrm)).format(n) for n in libraries]
        id_diffs_list= hapProbsAll(haplist=id_diffs, hap='id')

    except Exception as e:
        e.args = (error_in,)
        raise


    #print("id_diffs_list...............................................................", id_diffs_list)
    diffs_list = findDiff(inds=probs_list, posall=pos_list)
    id_diffs_list[id_diffs_list==-9]=np.nan

    dwins, twins = getWin(df=diffs_list, pos=pos_list,interval=interval)
    id_dwins, id_twins = getWin(df=id_diffs_list, pos=pos_list,interval=interval)

    with open("interval.txt", 'w') as f:
        print(interval,file=f)

    return dwins,twins,id_dwins,id_twins, np.ones(len(twins)) * chrm

def parallel_mergedchrm(libraries, totalch, interval, cores):
    print("Merging all chromosomes..")
    pool = mp.Pool(cores)
    res=[pool.apply_async(get_merged_chrm, args=((libraries, chrm, int(interval)))) for chrm in totalch]
    pool.close()
    pool.join()

    error_in="Can not process input data. Please check all input files. In particular, bedfile and bamfiles should have chromosomes listed as 1,2,... and the bedfile should be tab-separated. Please check the README.md for example input files."
    try:
        allf = [p.get() for p in res]
    except Exception as e:
        e.args = (error_in,)
        raise

    df0, tf0, id_df0, id_tf0 = allf[0][0], allf[0][1], allf[0][2], allf[0][3]
    chrmf=allf[0][4].tolist()
    for i in range(1,len(allf)):
        df0 = np.concatenate((df0, allf[i][0]), axis=0)
        tf0 = np.concatenate((tf0, allf[i][1]), axis=0)
        id_df0 = np.concatenate((id_df0, allf[i][2]), axis=0)
        id_tf0 = np.concatenate((id_tf0, allf[i][3]), axis=0)
        chrmf.extend(allf[i][4])

    np.savetxt(fname="chrm_list.csv", X=chrmf, delimiter=',')

    print("Finished merging all chromosomes...")
    return df0, tf0, id_df0, id_tf0, chrmf


def divergence_vcf(vcf_file, bedfile, target_ind, contam_ind):

    command = """bcftools view %s -R %s -s %s,%s -m2 -M2 -U -o intermediate1.vcf""" %(vcf_file,bedfile,target_ind,contam_ind)
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

    command = """bcftools view intermediate1.vcf -g ^miss -o intermediate2.vcf"""
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

    command = "bcftools query -f %CHROM\t%POS[\t%GT]\n intermediate2.vcf -o filtered_vcf.txt"
    process = subprocess.Popen(command.split(' '), stdout=subprocess.PIPE)
    output, error = process.communicate()

    divergence= nhFile(alldiff='filtered_vcf.txt')

    return divergence


def contamFile(infile, outfile, targets, idfile):
   cnf_u=pd.read_csv(infile, sep="\t",header=0,index_col=None)
   cnf_u=cnf_u.loc[cnf_u['name'].isin(targets),:]

   cnf = cnf_u.set_index('name')
   cnf=cnf.loc[targets]
   cnf['name']=cnf.index
   cnf.reset_index(drop=True, inplace=True)
   names=[]
   contam=[]
   for ch1 in range(len(cnf)-1):
       for ch2 in range(ch1+1,len(cnf)):
           names.append(cnf.loc[ch1, 'name'] + '_._' + cnf.loc[ch2, 'name'])
           ctotal=cnf.loc[ch1, 'contamination'] + cnf.loc[ch2, 'contamination']
           contam.append(ctotal)

   df = pd.DataFrame(
       {'name': names,
        'contam': contam
       })

   df_id = pd.DataFrame(
       {'name': cnf['name'],
        'contam': 2*cnf['contamination']
       })

   with pd.option_context('display.max_rows', len(df.index), 'display.max_columns', len(df.columns)):
   	    df.to_csv(outfile, sep=',')

   with pd.option_context('display.max_rows', len(df_id.index), 'display.max_columns', len(df_id.columns)):
   	    df_id.to_csv(idfile, sep=',')





def adjust_d(D_obs, c, p_c, p_e):
    x1 = p_e * (1-c)
    x2 = p_c * c
    return D_obs * x1 / (x1 + x2)

def adjust_s(S_obs, c, p_c, p_e):
    x1 = (1-p_e) * (1-c)
    x2 = (1-p_c) * c
    x=S_obs * x1 / (x1 + x2)
    x[np.isnan(x)]=0
    return x

def adjust_sd(D_obs, N_obs, *args, **kwargs):
    d_adj = adjust_d(D_obs, *args, **kwargs)
    s_adj = adjust_s(N_obs - D_obs, *args, **kwargs)
    return d_adj, d_adj+s_adj


def contamAll(diff, total, cfile, p_c):
    df=pd.read_csv(cfile, sep=",",header=0,index_col=0)
    cnt1=np.array(df['contam'])

    c=cnt1/100

    propmean=np.sum(diff, axis=0)/np.sum(total, axis=0)

    p_e = (propmean - c * p_c)/ (1-c)
    p_e[p_e>1]=1

    #Avoiding issues with low coverage libraries with contamination
    less0=np.where(p_e<0)[0]
    if len(less0)>0:
        print("Warning: In contamination correction step we detect some libraries with almost no data, or incorrect contamination levels. Data for comparison of these libraries is removed to avoid problems (see missing rows in hmm_parameters/p_all.csv).")
        diff[:,less0]=0
        total[:,less0]=0
        p_e[less0]=0
        c[less0]=0

    d_cor, n_cor = adjust_sd(D_obs=diff, N_obs=total, p_e=p_e, p_c=p_c, c=c)

    if np.all(np.where(np.isnan(np.sum(d_cor,0)))[0] == np.where(np.isnan(np.sum(n_cor,0)))[0]):
        d_cor[np.isnan(d_cor)] = 0
        n_cor[np.isnan(n_cor)] = 0

    if(np.sum(np.isnan(d_cor))==0):
        if(np.sum(np.isnan(n_cor))==0):
            return d_cor, n_cor
        else:
            print("Something is wrong in contamination correction, there are nan values.")


def get_length(x):
    try:
        return len(x)
    except TypeError:
        return 1

def data2p(diff_cor, total_cor, id_diff_cor, id_total_cor, libraries, listf, hmm_param, thresh, outdiff, outtotal, id_outdiff, id_outtotal, badwins):
    #print("diffs",diff_cor[[175,255]])
    #print("total",total_cor[total_cor<-1])
    print("Estimating p_0...")
    if get_length(badwins)==0:
        rem_wins=getHighDiv(alld=diff_cor, allt=total_cor)
    else:
        rem_wins=np.array(badwins).astype(int)
    np.savetxt(fname='filtered_windows.txt', X=rem_wins, delimiter=',')
    diff_cor[rem_wins,:] = 0
    total_cor[rem_wins,:] = 0

    id_diff_cor[rem_wins,:] = 0
    id_total_cor[rem_wins,:] = 0
    #print("diff",diff_cor[rem_wins,:])
    #print("total",total_cor[rem_wins,:])
    p1=getP(obsd=diff_cor, obst=total_cor, targets=listf,
        pfile=hmm_param+'p_0.txt',goodpairs='goodpairs.csv',
        allp=hmm_param+'p_all.csv', overf='overlap.csv',  thresh=thresh)

    getP(obsd=id_diff_cor, obst=id_total_cor, targets=libraries,
        pfile=hmm_param + 'p_2.txt',goodpairs='goodlibraries.csv',
        allp=hmm_param + 'identical_p_all.csv', overf='identical_overlap.csv',  thresh=thresh)

    np.savetxt(fname=outdiff, X=diff_cor, delimiter=',')
    np.savetxt(fname=outtotal, X=total_cor, delimiter=',')
    np.savetxt(fname=id_outdiff, X=id_diff_cor, delimiter=',')
    np.savetxt(fname=id_outtotal, X=id_total_cor, delimiter=',')

    return diff_cor, total_cor, id_diff_cor, id_total_cor, p1


def run_hmm(diff, total, chrm1, hbdf, likf, libraries, p1, cores):

    pool = mp.Pool(cores)
    [pool.apply_async(hmm, args=(diff, total, lib, np.array(chrm1), p1,
                    hbdf + 'pw_%s.csv' %(lib), likf + 'pw_%s.txt' %(lib), libraries)) for lib in libraries]
    pool.close()
    pool.join()


def pipeline1(targetsfile, bedfile, cores, rawbams, interval, splitbams, bedfiles, hapProbs, hmm_param, hbdf, likf, chrmf, sort1):

    libraries, listf = prep_function(targetsfile, splitbams, bedfiles, hapProbs, hmm_param, hbdf, likf)

    split_bed(bedfile=bedfile, totalch=chrmf, bedfiles=bedfiles)

    lib_chrm_all=create_lib_chrm(libraries=libraries, totalch=chrmf)
    if str(sort1)=='1':
        parallel_indexes(rawbams=rawbams, libraries=libraries, cores=cores)

    parallel_bamfilter(rawbams=rawbams, splitbams=splitbams, bedfiles=bedfiles, hapProbs=hapProbs, lib_chrm_all=lib_chrm_all, cores=cores, sort1=sort1)

    dwins,twins,id_dwins,id_twins, chrmlist = parallel_mergedchrm(libraries=libraries, totalch=chrmf, interval=interval, cores=cores)

    return libraries, listf, dwins, twins, id_dwins, id_twins, chrmlist
