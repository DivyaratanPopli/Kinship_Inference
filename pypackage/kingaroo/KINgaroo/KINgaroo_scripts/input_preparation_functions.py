#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 13:22:51 2019

@author: divyaratan_popli
"""


import numpy as np
import pandas as pd
import math
import scipy.special as scp
import pysam
import pybedtools

##########################################################start with simulation for relatedness distn



def pileup2alleles(bam, snp):
    """returns the number of ref and alt alleles at each position"""
    ref, alt = snp.name, snp.score
    ref_cnt, alt_cnt = 0, 0
    for pileupcolumn in bam.pileup(snp.chrom, snp.start, snp.end):
        if pileupcolumn.reference_pos == snp.start:     #find out the position in pilup that matches position in bed
            for base in pileupcolumn.get_query_sequences():     #get all bases present at that position
                if base.upper() == ref.upper(): ref_cnt +=1
                if base.upper() == alt.upper(): alt_cnt +=1
            return ref_cnt, alt_cnt
    return 0, 0


def makeHapProbs(inbam,inbai,inbed,outprob, outdiff_id):
    """takes bam, index, bed files to do random sampling and create pseudohaps"""
    diff_id=[]
    probsA=[]
    chrm=[]
    pos=[]

    bed = pybedtools.BedTool(inbed)
    with pysam.AlignmentFile(inbam) as bam:
        for i, snp in enumerate(bed):

            freq = pileup2alleles(bam, snp)
            refA=freq[0]
            altA=freq[1]
            nsites=refA+altA
            if nsites==0:
                probsA.append(-9)
                diff_id.append(-9)

            elif nsites>0 and nsites<2:

                p=altA/nsites
                probsA.append(p)

                diff_id.append(-9)


            elif nsites>=2:

                p=altA/nsites
                probsA.append(p)

                d0=scp.comb(altA,1) * scp.comb(refA,1)/scp.comb(nsites,2)
                diff_id.append(d0)

            chrm.append(snp.chrom)
            pos.append(snp.start)

    df=pd.DataFrame({
    "chrm":chrm,
    "pos":pos,
    "alt_prob":probsA
    })

    diffs=pd.DataFrame({
    "chrm":chrm,
    "pos":pos,
    "diff":diff_id
    })

    with pd.option_context('display.max_rows', len(df.index), 'display.max_columns', len(df.columns)):
            df.to_csv(outprob, sep=',')

    with pd.option_context('display.max_rows', len(diffs.index), 'display.max_columns', len(diffs.columns)):
            diffs.to_csv(outdiff_id, sep=',')



def hapProbsAll(haplist,hap):
    outarr=[]

    for fi in haplist:
        x=pd.read_csv(fi, sep=",", header=0,index_col=0)
        if hap=='id':
            x1=x['diff'].values
        elif hap=='noid':
            x1=x['alt_prob'].values

        outarr.append(x1)

    outarr=np.array(outarr).T


    if hap=='noid':

        pos=x['pos'].values
        chrm=x['chrm'].values

        return outarr, pos, chrm

    return outarr



def findDiff(inds,posall):
    """Find pairwise differences along genome for two specimen. In output file -9 is for missing data"""

    inds[inds==-9]=np.nan


    diff0=(np.vstack(inds[:,0]) * (1-inds)) + ((1-np.vstack(inds[:,0])) * inds)
    diffadd0=np.delete(diff0, np.s_[0:1], axis=1)
    for j in range(1,np.shape(inds)[1]):


        diffi=(np.vstack(inds[:,j]) * (1-inds)) + ((1-np.vstack(inds[:,j])) * inds)
        diffadd=np.delete(diffi, np.s_[0:j+1], axis=1)
        diffadd0=np.concatenate((diffadd0,diffadd),axis=1)

    return diffadd0


def getWin(df, pos, interval=1e7):

    length=pos[-1]

    bounds= np.arange(0,length,interval)
    snp_bin=np.digitize(pos,bounds,right=True)
    [uniq,ind]=np.unique(snp_bin,return_index=True)
    ind= np.append(ind,len(pos))
    ind1= ind[1:-1]

    dfsub=np.vsplit(df,ind1)
    ds=np.zeros([len(dfsub), np.shape(df)[1]])
    ts=np.zeros([len(dfsub), np.shape(df)[1]])
    for i in range(len(dfsub)):
        ds[i]=np.nansum(dfsub[i],0)
        ts[i]=np.sum(~np.isnan(dfsub[i]),0)

    winlen= len(ind)-1

    totalwin=math.ceil(length/interval)
    if winlen !=totalwin:
        missingwins=np.setdiff1d(list(range(1,totalwin+1)), uniq)
        ds=np.insert(ds,missingwins,np.zeros(np.shape(df)[1]),axis=0)
        ts=np.insert(ts,missingwins,np.zeros(np.shape(df)[1]),axis=0)


    return ds, ts


def nhFile(alldiff):

    df=pd.read_csv(alldiff, sep="\t",header=None,index_col=None, low_memory=False)

    df=df.loc[(df[2].apply(str)).isin(['0', '1']) & (df[3].apply(str)).isin(['0', '1']),:]

    nD,hD = np.array(df[2].apply(int)), np.array(df[3].apply(int))

    nA,hA = 1-nD, 1-hD

    diff=((nA*hD)+(nD*hA)) / ((nA+nD)*(hA+hD))
    nh_diff=np.mean(diff)

    return nh_diff


def getHighDiv(alld, allt):

    with np.errstate(divide='ignore', invalid='ignore'):
        prop=alld/allt
        avgprop= np.nansum(prop,axis=1) / np.sum(~np.isnan(prop),1)

    m=np.mean(avgprop[~np.isnan(avgprop)])
    sd=np.sqrt(np.var(avgprop[~np.isnan(avgprop)]))

    #print(m,sd)
    avgprop[np.isnan(avgprop)]=-9
    fres=np.where((avgprop>-9) & (avgprop>m+3*sd))[0]

    return fres


def getP(obsd, obst, targets, pfile, goodpairs, allp, overf, thresh=10):
    """Calculate the median proportion of differences in all pair of individuals. This will be an input for hmm"""
    med=[]
    names=[]

    with np.errstate(divide='ignore', invalid='ignore'):
        prop=np.sum(obsd, axis=0)/np.sum(obst, axis=0)

    goodlibs=np.sum(obst>0, axis=0)>thresh
    prop1=prop[goodlibs]

    names=np.array(targets)

    names1=names[goodlibs]


    med=np.median(prop1)


    over=np.sum(obst, axis=0)

    with open(pfile, 'w') as f:
        print(med,file=f)

    np.savetxt(fname=goodpairs, fmt='%s', X=names1, delimiter='/n')

    df = pd.DataFrame(
        {'pair': names,
        'prop': prop
        })

    with pd.option_context('display.max_rows', len(df.index), 'display.max_columns', len(df.columns)):
                df.to_csv(allp, sep=',')


    df_over=pd.DataFrame(
        {'pair': names,
        'overlap': over
        })

    with pd.option_context('display.max_rows', len(df_over.index), 'display.max_columns', len(df_over.columns)):
                df_over.to_csv(overf, sep=',')

    return med
