#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 13:22:51 2019

@author: divyaratan_popli
"""


import os
import msprime as msp
import numpy as np
import random
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib import figure
from scipy.stats import bernoulli
import math
import scipy.special as scp

#from astropy.stats import jackknife_resampling
#from astropy.stats import jackknife_stats
#from scipy.signal import savgol_filter
from ast import literal_eval
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


def hapProbsAll(haplist,hap,outf,posf,chrmf,psf):
    outarr=[]
    outps=[]

    for fi in haplist:
        x=pd.read_csv(fi, sep=",", header=0,index_col=0)
        if hap=='id':
            x1=x['diff'].values
        elif hap=='noid':
            x1=x['alt_prob'].values


            gudpos=np.where((x1<=1) & (x1>=0))
            x2=np.ones(np.shape(x1))*(-9)
            x2[gudpos]=np.random.binomial(1,x1[gudpos])
            outps.append(x2)

        outarr.append(x1)

    outarr=np.array(outarr).T
    np.savetxt(fname=outf, X=outarr, delimiter=',')

    if hap=='noid':

        pos=x['pos'].values
        chrm=x['chrm'].values
        np.savetxt(fname=posf, X=pos, delimiter=',')
        np.savetxt(fname=chrmf, X=chrm, delimiter=',')

        outps=np.array(outps).T
        np.savetxt(fname=psf, X=outps, delimiter=',')


def filterPseudohap(filist,posfile,outfilter,readfile,fil,fil_version):
    """Find polymorphic positions in all bam files listed in filist and store it in outfilter."""
    """Create per chromosome input files for READ in tped format"""

    poly=np.loadtxt(filist,dtype='float', delimiter = ",")
    pos=np.loadtxt(posfile,dtype='float', delimiter = ",")
    poly[poly==-9]=np.nan
    #ascpos=np.nansum(poly-np.vstack(poly[:,0]), axis=1)!=0
    if int(fil)==1:
        ascpos=[]
        for i in range(len(poly)):
            ascpos.append(len(np.unique(poly[i,~np.isnan(poly[i,:])]))>1)

        keep=pos[ascpos]


    elif int(fil)==0:
        keep=pos


    np.savetxt(fname=outfilter, X=keep, delimiter=",")

    if fil_version=='hap':
        read=np.repeat(poly, 2, axis=1)

        if int(fil)==1:
            readpos=read[ascpos]
        elif int(fil)==0:
            readpos=read

        np.savetxt(fname=readfile, X=readpos, delimiter=",")


def READInput(filist, posfile, read, tfam, libraries):
    """Merge per chromosome files from filterPseudohap into both filtered/non filtered tped files for READ, also output tfam file for READ."""



    fi0=np.loadtxt(filist[0],dtype='float', delimiter=",")
    po0=np.loadtxt(posfile[0],dtype='float', delimiter=",")

    chrm=[1]*len(po0)

    for i in range(1,len(filist)):

        fi=np.loadtxt(filist[i],dtype='float', delimiter=",")
        fi0=np.concatenate((fi0,fi),axis=0)

        po=np.loadtxt(posfile[i],dtype='float', delimiter=",")
        po0=np.concatenate((po0,po),axis=0)

        chrm.extend([i+1]*len(po))

    df=pd.DataFrame(columns=list(range(np.shape(fi0)[1])), index=list(range(np.shape(fi0)[0])), data=fi0)
    df[df==0]='A'
    df[df==1]='G'
    df=df.replace(np.nan,0)


    df.insert(loc=0, column='chrm', value=chrm)
    df.insert(loc=1, column='pos', value=po0)
    df.insert(loc=1, column='gndlist', value=0)
    df.insert(loc=1, column='name', value=df['chrm'].astype(str)+'_'+df['pos'].astype(str))
    df['pos']=df['pos'].astype(int)


    with pd.option_context('display.max_rows', len(df.index), 'display.max_columns', len(df.columns)):
            df.to_csv(read, sep='\t', header=None, index=False)


    #making the tfam files
    libs=['_'+str(i)+'_' for i in libraries]
    tf=pd.DataFrame()
    tf[0]=libs
    tf[1]=libs
    tf[2]=0
    tf[3]=0
    tf[4]=0
    tf[5]=0

    with pd.option_context('display.max_rows', len(tf.index), 'display.max_columns', len(tf.columns)):
            tf.to_csv(tfam, sep='\t', header=None, index=False)





def findDiff(indfile,infilter,posall,difffile,fil):
    """Find pairwise differences along genome for two specimen. In output file 9 is for missing data"""
    """Fil=1 is for using only polymorphic sites"""

    inds=np.loadtxt(indfile, dtype='float', delimiter=",")

    if int(fil)==1:
        pos_fil=np.loadtxt(infilter, dtype='float', delimiter=",")
        pos_all=np.loadtxt(posall, dtype='float', delimiter=",")

        filIdx = np.isin(pos_all,pos_fil)

        inds1=inds[filIdx]

    elif int(fil)==0:
        inds1=inds.copy()

    inds1[inds1==-9]=np.nan


    diff0=(np.vstack(inds1[:,0]) * (1-inds1)) + ((1-np.vstack(inds1[:,0])) * inds1)
    diffadd0=np.delete(diff0, np.s_[0:1], axis=1)
    for j in range(1,np.shape(inds1)[1]):


        diffi=(np.vstack(inds1[:,j]) * (1-inds1)) + ((1-np.vstack(inds1[:,j])) * inds1)
        diffadd=np.delete(diffi, np.s_[0:j+1], axis=1)
        diffadd0=np.concatenate((diffadd0,diffadd),axis=1)

    np.savetxt(fname=difffile, X=diffadd0, delimiter=",")




def identicalDiff(indfile, infilter, posall, difffile, fil):
    inds=np.loadtxt(indfile, dtype='float', delimiter=",")

    if int(fil)==1:
        pos_fil=np.loadtxt(infilter, dtype='float', delimiter=",")
        pos_all=np.loadtxt(posall, dtype='float', delimiter=",")

        filIdx = np.isin(pos_all,pos_fil)

        inds1=inds[filIdx]
    elif int(fil)==0:
        inds1=inds.copy()

    inds1[inds1==-9]=np.nan
    np.savetxt(fname=difffile, X=inds1, delimiter=",")



def getWin(diff, posf, dfile, tfile, interval):
    df=np.loadtxt(diff, dtype='float', delimiter=',')
    pos=np.loadtxt(posf, dtype='float', delimiter=',')

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


    np.savetxt(fname=dfile, X=ds, delimiter=",")
    np.savetxt(fname=tfile, X=ts, delimiter=",")


def mergeChrm(dinfiles, tinfiles, dfile, tfile, chfile):
    """Merge per chromosome window files for a pair of individual"""
    #merging the data files


    df0=np.loadtxt(dinfiles[0], dtype='float', delimiter=',')
    tf0=np.loadtxt(tinfiles[0], dtype='float', delimiter=',')
    chrm0=np.ones(len(tf0))
    for fi in range(1,len(dinfiles)):

        df=np.loadtxt(dinfiles[fi], dtype='float', delimiter=',')
        tf=np.loadtxt(tinfiles[fi], dtype='float', delimiter=',')
        chrm=np.ones(len(tf)) * (fi+1)

        df0=np.concatenate((df0,df),axis=0)
        tf0=np.concatenate((tf0,tf),axis=0)
        chrm0=np.concatenate((chrm0,chrm),axis=0)

    np.savetxt(fname=dfile, X=df0, delimiter=",")
    np.savetxt(fname=tfile, X=tf0, delimiter=",")
    np.savetxt(fname=chfile, X=chrm0, delimiter=",")

def contamFile(infile, outfile, targets, idfile):

    cnf=pd.read_csv(infile, sep="\t",header=0,index_col=None)
    cnf=cnf.loc[cnf['name'].isin(targets),:]
    cnf.reset_index(drop=True, inplace=True)
    names=[]
    contam=[]
    for ch1 in range(len(cnf)-1):
        for ch2 in range(ch1+1,len(cnf)):
            names.append(cnf.loc[ch1, 'name'] + '_' + cnf.loc[ch2, 'name'])
            ctotal=cnf.loc[ch1, 'contamination (%)'] + cnf.loc[ch2, 'contamination (%)']
            contam.append(ctotal)

    df = pd.DataFrame(
        {'name': names,
         'contam': contam
        })

    df_id = pd.DataFrame(
        {'name': cnf['name'],
         'contam': 2*cnf['contamination (%)']
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


def contamAll(dfile, tfile, cfile, Dnhfile, difffile, totalfile, iscnt):

    diff=np.loadtxt(dfile, delimiter=',', dtype='float')
    total=np.loadtxt(tfile, delimiter=',', dtype='float')

    if str(iscnt)=='1':
        p_c=Dnhfile


        df=pd.read_csv(cfile, sep=",",header=0,index_col=0)
        cnt1=np.array(df['contam'])

        c=cnt1/100

        propmean=np.sum(diff, axis=0)/np.sum(total, axis=0)


        p_e = (propmean - c * p_c)/ (1-c)
        p_e[p_e>1]=1

        d_cor, n_cor = adjust_sd(D_obs=diff, N_obs=total, p_e=p_e, p_c=p_c, c=c)
        #print(d_cor[np.isnan(d_cor)])

        if np.all(np.where(np.isnan(np.sum(d_cor,0)))[0] == np.where(np.isnan(np.sum(n_cor,0)))[0]):
            d_cor[np.isnan(d_cor)] = 0
            n_cor[np.isnan(n_cor)] = 0


        if(np.sum(np.isnan(d_cor))==0):
            if(np.sum(np.isnan(n_cor))==0):

                np.savetxt(fname=difffile, X=d_cor, delimiter=",")
                np.savetxt(fname=totalfile, X=n_cor, delimiter=",")

    elif str(iscnt)=='0':
        np.savetxt(fname=difffile, X=diff, delimiter=",")
        np.savetxt(fname=totalfile, X=total, delimiter=",")

def getHighDiv(diff, total, fout):

    alld= np.loadtxt(diff, delimiter=',', dtype='float')
    allt= np.loadtxt(total, delimiter=',', dtype='float')


    prop=alld/allt
    avgprop= np.nansum(prop,axis=1) / np.sum(~np.isnan(prop),1)

    m=np.mean(avgprop[~np.isnan(avgprop)])
    sd=np.sqrt(np.var(avgprop[~np.isnan(avgprop)]))

    #print(m,sd)
    avgprop[np.isnan(avgprop)]=-9
    fres=np.where((avgprop>-9) & (avgprop>m+3*sd))[0]
    np.savetxt(fout, fres, delimiter=',')


def remHighDiv(diffname, totalname, winfile, outd, outt):
    winsHD=np.loadtxt(winfile, ndmin=1).astype(int)

    df= np.loadtxt(diffname, delimiter=',', dtype='float')
    dt= np.loadtxt(totalname, delimiter=',', dtype='float')
    df[winsHD,:] = 0
    dt[winsHD,:] = 0

    np.savetxt(outd, df, delimiter=',')
    np.savetxt(outt, dt, delimiter=',')


def getP(dfile, tfile, targets, pfile, goodpairs, allp, overf):
    """Calculate the median proportion of differences in all pair of individuals. This will be an input for hmm"""
    med=[]
    names=[]



    obsd= np.loadtxt(dfile, delimiter=',', dtype='float')
    obst= np.loadtxt(tfile, delimiter=',', dtype='float')

    prop=np.sum(obsd, axis=0)/np.sum(obst, axis=0)
    #chrm= np.loadtxt(chfile, delimiter=',', dtype='float')
    #calculating p1
    goodlibs=np.sum(obst>0, axis=0)>10

    obsd1=obsd[:,goodlibs]
    obst1=obst[:,goodlibs]
    prop1=prop[goodlibs]

    names=np.array(targets)

    names1=names[goodlibs]


    med=np.median(prop1)


    over=np.sum(obst, axis=0)

    with open(pfile, 'w') as f:
        print(med,file=f)
    #np.savetxt(fname=allp, X=med, delimiter=',')
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


def getOverlap(filist, overfile):

    cols=["pair", "overlap"]
    df=pd.DataFrame(columns=cols)
    for i,fi in enumerate(filist):
        pair=fi.split('pw_')[1].split('.cs')[0]
        data=pd.read_csv(fi, sep=',',index_col=0,header=0)
        df.loc[i,"overlap"]=np.sum(data['count'])
        df.loc[i,"pair"]=pair

    with pd.option_context('display.max_rows', len(df.index), 'display.max_columns', len(df.columns)):
            df.to_csv(overfile, sep=',')



def inbredwin_transform(inputall,alls,targets,fi,outlist,totalch,interval,fil):

    tb_all=pd.read_csv(inputall,sep=',',header=0,index_col=0)
    unique_id=np.unique(tb_all['id'])
    l=len(fi)
    allids=np.array(list(range(1,l+1)))
    noinbred=np.setdiff1d(allids,unique_id,assume_unique=False)


    totalch=np.array(totalch)
    print(totalch[0])
    for idn in unique_id:
        tb=tb_all.loc[tb_all['id']==idn,:]


        allsites=pd.read_csv(alls,sep='\t',header=None,index_col=False, low_memory=False)
        allsites=allsites.loc[allsites[0].isin(totalch.astype(str)),:]
        allsites[0]=allsites[0].astype(int)

        ind=allsites[0].diff()[allsites[0].diff() != 0].index.values - 1
        ind=ind[1:]
        ind=np.append(ind,len(allsites[0])-1)


        wint=[]

        for ch in totalch:

            length=allsites.loc[ind[ch - 1],1]

            inbred=tb.loc[tb['chrom']==ch,['start_pos','end_pos']]
            for i in range(np.shape(inbred)[0]):

                pos=inbred.iloc[i,:]

                bounds= np.arange(0,length,interval)
                snp_bin=np.digitize(pos,bounds,right=True)
                print(snp_bin)
                if snp_bin[0]==snp_bin[1]:
                    bins=snp_bin[0]
                else :
                    bins=list(range(snp_bin[0],snp_bin[1]+1))
                wint.append((ch,bins))



        pair=targets[idn-1]


        obs=pd.read_csv(fi[idn-1],sep=',',header=0,index_col=0)
        data=obs[['dis','count']]


        inbredwin=[]


        for i in range(len(wint)):
            winlist=[]

            if type(wint[i][1]) is list:
                winlist.extend(wint[i][1])
            else:
                winlist.append(wint[i][1])

            in_ch=wint[i][0]
            for j in range(len(winlist)):
                in_win=winlist[j]


                zero=obs.loc[obs["chrom"]==in_ch,:].index.values[0]-1
                inbredwin.append(zero+in_win)

        hbdstate=np.ones(len(data))
        hbdstate[inbredwin]=0
        print(hbdstate)
        np.savetxt(fname=outlist[idn-1], X=hbdstate)


    for noin in noinbred:
        pair=targets[noin-1]

        obs=pd.read_csv(fi[noin-1],sep=',',header=0,index_col=0)
        data=obs[['dis','count']]
        hbdstate=np.ones(len(data))
        np.savetxt(fname=outlist[noin-1], X=hbdstate)
