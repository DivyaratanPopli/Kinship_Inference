#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 21:36:58 2021

@author: divyaratan_popli
"""


#import scipy
import numpy as np
import random
import pandas as pd
#from scipy.stats import binom
#from math import isclose
#from scipy.special import gammaln
#from scipy import optimize
#import math
#import matplotlib.pyplot as plt
#from scipy.special import loggamma
#from scipy.special import beta as Betaf
#from scipy.optimize import minimize_scalar
#import pylab
#from numpy.linalg import matrix_power
#from scipy.special import logsumexp
#import numba


def getpc_sib(filist, pairs, truerel,kinrel_nf,pclist, siblist, grlist,halfsib,avulist,thirdlist,fourthlist,fifthlist,idlist,runtable,ctw):

    cols=["pair", "relatedness","second_guess", "loglik_ratio"]
    df=pd.DataFrame(columns=cols)
    runtable_cut=pd.DataFrame(0,index=truerel, columns=kinrel_nf)

    for i,fi in enumerate(filist):
        pair=pairs[i]
        lik=np.loadtxt(fi,dtype='float', delimiter = ",")
        data=pd.DataFrame(lik.reshape(-1,len(truerel)), columns=truerel)

        b=data.sort_values(by=0, ascending=False, axis=1)

        print(b)


        df.loc[i,"relatedness"]=b.columns[0]
        df.loc[i,"pair"]=pair

        df.loc[i,"second_guess"]=b.columns[1]
        df.loc[i,"loglik_ratio"]=(b[b.columns[0]] - b[b.columns[1]]).item()

        l=df.loc[i,"pair"]

        if l in grlist:
            label='gr'
        elif l in pclist:
            label='pc'
        elif l in siblist:
            label='sib'
        elif l in idlist:
            label='id'
        elif l in thirdlist :
            label='deg3'
        elif l in fourthlist :
            label='deg4'
        elif l in halfsib :
            label='hsib'
        elif l in avulist :
            label='avu'
        elif l in fifthlist :
            label='deg5'
        else:
            label='un'


        if df.loc[i,'loglik_ratio'] >= float(ctw):
            runtable_cut.loc[label,df.loc[i,"relatedness"]]=runtable_cut.loc[label,df.loc[i,"relatedness"]]+1
        else:
            runtable_cut.loc[label,'NF']=runtable_cut.loc[label,'NF']+1

    with pd.option_context('display.max_rows', len(runtable_cut.index), 'display.max_columns', len(runtable_cut.columns)):
            runtable_cut.to_csv(runtable, sep=',')




def bigTable(tables, alltable):
    filenames=tables
    lik_colind=pd.read_csv(filenames[0], sep=",", header=0,index_col=0)
    c=lik_colind.columns.values
    r=lik_colind.index.values
    bigtable=pd.DataFrame(0,index=r, columns=c)
    for t in tables:
        lik=pd.read_csv(t, sep=",", header=0,index_col=0)
        bigtable=bigtable+lik
    bigtable["sum"] = bigtable.sum(axis=1)
    bigtable = bigtable.loc[:,c].div(bigtable["sum"], axis=0)

    with pd.option_context('display.max_rows', len(bigtable.index), 'display.max_columns', len(bigtable.columns)):
        bigtable.to_csv(alltable, sep=',')
