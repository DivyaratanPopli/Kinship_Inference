#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 21:36:58 2021

@author: divyaratan_popli
"""


import scipy
import numpy as np
import random
import pandas as pd
from scipy.stats import binom
from math import isclose
from scipy.special import gammaln
from scipy import optimize
import math
import matplotlib.pyplot as plt
from scipy.special import loggamma
from scipy.special import beta as Betaf
from scipy.optimize import minimize_scalar
import pylab
from numpy.linalg import matrix_power
from scipy.special import logsumexp
import numba

#print(scipy.__version__)



"""
#to stop the program if there is a warning
import warnings

with warnings.catch_warnings():
    warnings.simplefilter('error')
    for chrm in range(chrm_no):

            for t in range(pos[chrm],pos[chrm+1]-1):


                for i in range(M):
                    print.info(chrm,t,i)
                    numerator_l = np.log(alpha[t, i]) + np.log(A[i, :]) + np.log(B[:, t + 1].T) + np.log(beta[t + 1, :].T) - np.log(scale[t+1])  #it is the probability that hmm path has particular arc at time t
                    xi_l[i, :, t] = numerator_l


"""
@numba.vectorize
def Betal(x, y):
    # numba doesn't understand scipy.special.gammaln, but it understands math.lgamma,
    # which is the same, except that it doesn't understand arrays.
    return math.lgamma(x) + math.lgamma(y) - math.lgamma(x+y)

@numba.vectorize
def fact(x):
    return math.lgamma(x + 1)

@numba.njit
def forward(A,B,pi,data):
    N = A.shape[0]
    T = data.shape[0]
    alpha = np.zeros((T,N))
    alpha_p = np.zeros((T,N))
    scale = np.zeros(T)

    alpha[0,:] = pi * B[:,0]
    alpha_p[0,:] = alpha[0,:] / alpha[0,:].sum()
    scale[0] = alpha[0,:].sum()

    for t in range(1,T):
        alpha[t,:]= (alpha_p[t-1,:] @ A) * B[:,t]
        alpha_p[t,:]= alpha[t,:] / alpha[t,:].sum()
        scale[t] = alpha[t,:].sum()

    return alpha_p, scale

#@numba.njit
def backward(A, B, data, scale):

    N= np.shape(A)[0]
    T= np.shape(data)[0]
    beta=np.zeros((T,N))
    beta_p=np.zeros((T,N))

    beta[-1,:]=np.ones(N)
    beta_p[-1,:]=np.ones(N)

    for t in reversed(range(T-1)):
        for n in range(N):
            beta[t,n]= sum( A[n,:] * beta_p[t+1,:] * B[:,t+1])
            beta_p[t,n]= beta[t,n]/scale[t+1]

    return beta_p



@numba.njit
def objective(a, meanp, k, data):
    cost=0

    b=a * (1-meanp)/meanp

    for w in range(len(data)):

        ll=(fact(data[w,1]) - fact(data[w,0]) - fact(data[w,1]-data[w,0]) + \
             Betal(data[w,0]+a, data[w,1]-data[w,0]+b) - Betal(a,b))


        Wll= ll * k[w]

        cost= cost + Wll

    return (-1)*cost



def makeB(data,xin,M,T,inbr):

    B = np.zeros((inbr,M,T))

    for st in range(inbr):
        for re in range(M):
            a=xin[st,re,0]
            b=xin[st,re,1]
            B[st,re,:]=np.exp(fact(data[:,1]) - fact(data[:,0]) - fact(data[:,1]-data[:,0]) + Betal(data[:,0]+a, data[:,1]-data[:,0]+b) - Betal(a,b))

    return B

def getEmissions(B, hbd):
    B2=np.zeros([np.shape(B)[1],np.shape(B)[2]])
    scale1=np.zeros(np.shape(B)[2])

    for i in range(np.shape(B)[1]):
        B2[i,:]=(B[0,i,:] * hbd[:,0]) + (B[1,i,:] * hbd[:,1]) + (B[2,i,:] * hbd[:,2])

    scale1=B2.max(axis=0)
    B2scaled=B2/scale1

    return B2scaled, np.sum(np.log(scale1))


def expectation(data, hbd, A, B, pos, chrm_no, pi):

    M = A.shape[0]
    T = len(data)
    alpha=np.zeros((T,M))
    beta=np.zeros((T,M))
    scale=np.zeros(T)

    B2d,log_bscale1 = getEmissions(B, hbd)


    for chrm in range(chrm_no):
        data_chr=data[pos[chrm]:pos[chrm+1],:]

        alpha[pos[chrm]:pos[chrm+1],:],scale[pos[chrm]:pos[chrm+1]] = forward(A,B2d[:,pos[chrm]:pos[chrm+1]],pi,data_chr)
        beta[pos[chrm]:pos[chrm+1],:] = backward(A,B2d[:,pos[chrm]:pos[chrm+1]],data_chr,scale[pos[chrm]:pos[chrm+1]])

    post= alpha*beta

    li=np.sum(np.log(scale))+ log_bscale1


    assert isclose(sum(abs(np.sum(post,axis=1) - np.ones(np.shape(B2d)[1]))), 0, abs_tol=1e-7), "sum of alpha and beta is not 1"
    assert isclose(sum(abs(np.sum(A,axis=1) - np.ones(np.shape(A)[0]))), 0, abs_tol=1e-7), "sum of transition matrix columns is not 1"

    pi_new=matrix_power(A,10000)[0]
    #pi_new=np.mean(post[pos[:-1],:],0)

    return post.T, li, pi_new


def makexin(x):
    xin=np.zeros([3,3,2])
    xin[0,:,:]=np.reshape(x[0:6],(3,2))
    xin[1,:,:]=[x[4:6],x[2:4],x[6:8]]
    xin[2,:,:]=[x[6:8],x[2:4],x[6:8]]
    #print("xin=",xin)
    return xin


def emission(gamma, data, hbd, p1avg, A, inbr, x0, bnds):

    M = A.shape[0]
    T = len(data)

    k_pc= gamma[0,:] * hbd[:,0]
    k_un= gamma[1,:]
    k_id= (gamma[2,:] * hbd[:,0]) + (gamma[0,:] * hbd[:,1])
    k_inb= (gamma[2,:] * hbd[:,1]) + (gamma[2,:] * hbd[:,2]) + (gamma[0,:] * hbd[:,2])

    [mean_pc,mean_un,mean_id,mean_inb]=[p1avg[0,0],p1avg[0,1],p1avg[0,2],p1avg[2,2]]
    [a_pc, a_un, a_id, a_inb]=[x0[0],x0[2],x0[4],x0[6]]





    solution_pc = minimize_scalar(objective,a_pc,args=(mean_pc,k_pc,data),method='Bounded',bounds=[bnds[0],10000],options={'disp':1})
    a_pc=solution_pc.x
    b_pc= a_pc * (1-mean_pc)/mean_pc

    solution_un = minimize_scalar(objective,a_un,args=(mean_un,k_un,data),method='Bounded',bounds=[bnds[1],10000],options={'disp':1})
    a_un=solution_un.x
    b_un= a_un * (1-mean_un)/mean_un

    solution_id = minimize_scalar(objective,a_id,args=(mean_id,k_id,data),method='Bounded',bounds=[bnds[2],10000],options={'disp':1})
    a_id=solution_id.x
    b_id= a_id * (1-mean_id)/mean_id

    solution_inb = minimize_scalar(objective,a_inb,args=(mean_inb,k_inb,data),method='Bounded',bounds=[bnds[3],10000],options={'disp':1})
    a_inb=solution_inb.x
    b_inb= a_inb * (1-mean_inb)/mean_inb



    x=[a_pc,b_pc,a_un,b_un,a_id,b_id,a_inb,b_inb]

    xin=makexin(x)


    up_p=np.zeros([np.shape(xin)[0], np.shape(xin)[1]])
    for pi in range(np.shape(xin)[0]):
        for ps in  range(np.shape(xin)[1]):
            up_p[pi,ps]=xin[pi,ps,0] / (xin[pi,ps,0] + xin[pi,ps,1])


    Bu=makeB(data,xin,M,T,inbr)

    return Bu, up_p, xin


def bnd_calc(mean_p, dist, propvar):

    r=(1-mean_p)/mean_p
    t=propvar*dist*dist

    low_bnd=-1 * (t + (t*r*r) + r*((2*t)-1))/(t+ (3*t*r) + (3*r*r*t) + (r*r*r*t))
    return low_bnd

def baum_welch(data, hbd, A, B, pos, p1avg, inbr, x0, max_iter=1000):
    #print(data)
    chrm_no=len(pos)-1
    n_iter=0
    diff=1
    last_lik=0


    #print("initial beta parameters are:",x0)

    [mean_pc,mean_un,mean_id,mean_inb]=[p1avg[0,0],p1avg[0,1],p1avg[0,2],p1avg[2,2]]

    di=(mean_un-mean_pc)

    bnd_pc= bnd_calc(mean_p=mean_pc, dist=di, propvar=3/4)
    bnd_un= bnd_calc(mean_p=mean_un, dist=di, propvar=1)
    bnd_id= bnd_calc(mean_p=mean_id, dist=di, propvar=1/2)
    bnd_inb= bnd_calc(mean_p=mean_inb, dist=di, propvar=1)

    pi=np.array([1/3,1/3,1/3])

    while diff > 0.000001 and n_iter<max_iter:


        #print(pi)
        [gamma, li, pi]=expectation(data=data, hbd=hbd, A=A, B=B, pos=pos, chrm_no=chrm_no, pi=pi)

        [B, up_p, x1]=emission(gamma=gamma, data=data, hbd=hbd, p1avg=p1avg, A=A, inbr=inbr, x0=x0, bnds=[bnd_pc, bnd_un, bnd_id, bnd_inb])
        x0=np.reshape(x1,(1,18))[0,[0,1,2,3,4,5,10,11]]

        if n_iter !=0:
            diff= li - last_lik
            last_lik=li
        elif n_iter==0:
            diff= last_lik - li
            last_lik=li
        n_iter=n_iter+1
        #print("diff=",diff)
        #print("up_p=",up_p)

        #if diff<0:
            #print("likelihood is going negative. abort")

        #print("lik=",li)

    B2final,log_bfinal1=getEmissions(B, hbd)
    return gamma, A , B2final, up_p, li, pi


def viterbi(data, A, B, pi):

    T = data.shape[0]
    M = A.shape[0]

    omega = np.zeros((T, M))
    omega[0, :] = np.log(pi * B[:, 0])

    prev = np.zeros((T - 1, M))

    for t in range(1, T):
        for j in range(M):
            # Same as Forward Probability
            probability = omega[t - 1] + np.log(A[:, j]) + np.log(B[j,t])

            # This is our most probable state given previous state at time t (1)
            prev[t - 1, j] = np.argmax(probability)

            # This is the probability of the most probable state (2)
            omega[t, j] = np.max(probability)

    # Path Array
    S = np.zeros(T)

    # Find the most probable last hidden state
    last_state = np.argmax(omega[T - 1, :])


    S[0] = last_state

    backtrack_index = 1
    for i in range(T - 2, -1, -1):
        S[backtrack_index] = prev[i, int(last_state)]
        last_state = prev[i, int(last_state)]
        backtrack_index += 1

    # Flip the path array since we were backtracking
    S = np.flip(S, axis=0)
    S[S==2]=1/2
    S[S==0]=3/4
    return S

def hmm(difffile, targets, totalfile, listind, listf, pfile, Afiles, resfiles, gammafiles, likfile, pairf1, pairf2, instates, p1thresh, thresh, rels):

    inbr_states=len(instates)
    ind=np.where(listf==listind)[0][0]

    i1=np.where(targets.astype(str)==listind.split('_')[0])[0][0]
    i2=np.where(targets.astype(str)==listind.split('_')[1])[0][0]

    diff= np.loadtxt(difffile,dtype='float', delimiter = ",")[:,ind]
    total= np.loadtxt(totalfile,dtype='float', delimiter = ",")[:,ind]

    data=np.vstack([diff,total]).T
    if np.sum(data[:,1]>0)>10:
        pth=np.loadtxt(p1thresh,dtype='float', delimiter = ",")[:,[i1,i2]]


        pair1=pd.read_csv(pairf1, sep=",", header=0,index_col=0)
        pair2=pd.read_csv(pairf2, sep=",", header=0,index_col=0)

        pair1['infocount']=pth[:,0]
        pair2['infocount']=pth[:,1]
        ####correct for numeric errors from hbd_hmm
        pair1.loc[(pair1['g_noin']!=9) & (pair1['g_noin']>1),'g_noin']=1
        pair2.loc[(pair2['g_noin']!=9) & (pair2['g_noin']>1),'g_noin']=1

        pair1.loc[pair1['g_noin']<0,'g_noin']=0
        pair2.loc[pair2['g_noin']<0,'g_noin']=0

    #    pair1['ninb']=pair1['dis']/pair1['count']     #1 means no hbd
    #    pair1['inb']=1-pair1['ninb']
    #    pair2['ninb']=pair2['dis']/pair2['count']     #1 means no hbd
    #    pair2['inb']=1-pair2['ninb']

        pair1.loc[pair1['infocount']<thresh, 'g_noin']=9
        pair2.loc[pair2['infocount']<thresh, 'g_noin']=9

    #    pair1.loc[pair1['g_noin']>0.5,'g_noin']=1
    #    pair2.loc[pair2['g_noin']>0.5,'g_noin']=1

        pair1['ninb']=pair1['g_noin']     #1 means no hbd
        pair1['inb']=1-pair1['g_noin']
        pair2['ninb']=pair2['g_noin']     #1 means no hbd
        pair2['inb']=1-pair2['g_noin']


        badwin1=list(pair1.index[pair1['ninb']==9])
        badwin2=list(pair2.index[pair2['ninb']==9])

        pair1.loc[badwin1,'inb']=0
        pair1.loc[badwin1,'ninb']=1
        pair2.loc[badwin2,'inb']=0
        pair2.loc[badwin2,'ninb']=1



        hbd=np.zeros([len(total),3])
        hbd[:,0]=pair1['ninb']*pair2['ninb']
        hbd[:,1]=(pair1['ninb']*pair2['inb']) + (pair2['ninb']*pair1['inb'])
        hbd[:,2]=pair1['inb']*pair2['inb']

        gudwin=total>0
        data=data[gudwin,:]
        hbd=hbd[gudwin,:]

        chrmlist=np.asarray(pair1['chrom'])
        chrmlist=chrmlist[gudwin]

        win=np.sum(gudwin)

        with open(pfile,"r") as f:
            p_1=float(f.read())

        p_12=p_1/2
        inerror=p_12/2
        initial_p=np.array([[(p_1+p_12)/2,p_1,p_12],
                            [p_12,p_1,inerror],
                            [inerror,p_1, inerror]])


        pos=np.where(chrmlist[:-1] != chrmlist[1:])[0]+1
        pos=np.append(pos,np.shape(total)[0])
        pos=np.insert(pos, 0, 0, axis=0)



        #pi= np.array([1/3,1/3,1/3])

        #A=A1
        likall=[]
        for rel_cnt in range(len(rels)):

            Afile=Afiles[rel_cnt]
            A=np.array(pd.read_csv(Afile,sep=',',header=None,index_col=False))
            pi=matrix_power(A,10000)[0]
            Bu = np.zeros((inbr_states,np.shape(A)[0],win)) #Emission probability
            #for st in instates:
             #   for re in range(np.shape(A)[0]):
             #       Bu[st,re,:]=binom.pmf(data[:,0], data[:,1], initial_p[st,re])

            [b0, b1, b2, b3]=[1000,1000,1000,1000]
            [a0,a1,a2,a3]=[b0*initial_p[0,0]/(1-initial_p[0,0]), b1*initial_p[0,1]/(1-initial_p[0,1]), b2*initial_p[0,2]/(1-initial_p[0,2]), b3*initial_p[2,2]/(1-initial_p[2,2])]
            x0 = [a0,b0,a1,b1,a2,b2,a3,b3]

            xin0=makexin(x=x0)
            Bu=makeB(data=data,xin=xin0, M=A.shape[0], T=len(data), inbr=inbr_states)



            #print(A)
            #[B, log_bscale]= scaleB(Bu, instates)
            B=Bu

            np.argwhere(np.isnan(B))



            gamma,A,B,up_p,lik, pi= baum_welch(data=data, hbd=hbd, A=A, B=B, pos=pos, p1avg=initial_p, inbr=inbr_states, x0=x0)
            res=viterbi(data, A, B, pi)
            likall.append(lik)
            np.savetxt(fname=resfiles[rel_cnt], X=res,delimiter=',')
            np.savetxt(fname=gammafiles[rel_cnt], X=gamma.T,delimiter=',')

    else:
        likall=list([float('-inf')]*len(rels))
        for rel_cnt in range(len(rels)):
            res=np.ones(len(data))*-9
            gamma=np.ones([len(data),3])*-9

            np.savetxt(fname=resfiles[rel_cnt], X=res, delimiter=',')
            np.savetxt(fname=gammafiles[rel_cnt], X=gamma.T, delimiter=',')

    np.savetxt(fname=likfile, X=likall,delimiter=',')



"""

p1="Chagyrskaya1141"
p2="Chagyrskaya13"
pair=p1+"_"+p2
datafolder="/mnt/diversity/divyaratan_popli/100arc/inbreeding/laurits_writeup/Ch9_rem/neavar/noCh12_bayesContam/"
hbdfolder="/mnt/diversity/divyaratan_popli/100arc/inbreeding/laurits_writeup/Ch9_rem/neavar/noCh12_26mar_fil112_114/"
rel='sib'
fil='0'
pfile=datafolder+"hmm_parameters_fil%s/p1_file" %(fil)
phalf=datafolder+"hmm_parameters_fil%s/phalf_file" %(fil)
datafile=datafolder+"mergedwin_remHighdiv_fil%s/pw_%s.csv" %(fil,pair)
Afile="/mnt/diversity/divyaratan_popli/100arc/inbreeding/simulation_markov_inbreeding3/new_transitions1e6/transition_matrix_%s.csv" %(rel)
pairf1=hbdfolder+"hbd_win_fil%s/pw_%s.csv" %(fil,p1)
pairf2=hbdfolder+"hbd_win_fil%s/pw_%s.csv" %(fil,p2)

p1thresh=datafolder+"identicalmergedwin_fil1/pw_%s.csv" %(p1)
p2thresh=datafolder+"identicalmergedwin_fil1/pw_%s.csv" %(p2)
thresh=20
instates=[0,1,2]
inbr_states=3
inerror=0.001
upA=rel

plt.plot(data[:,0]/data[:,1])
plt.hlines(y=p_1,xmin=0,xmax=250)
plt.hlines(y=p_12,xmin=0,xmax=250)


d='0.15'
inb="1"
p1="4"
p2="11"
pair=p1+"_"+p2
run=15
datafolder="/mnt/diversity/divyaratan_popli/100arc/inbreeding/simulation_markov_inbreeding3/inbred%s/run%s/downsample%s/" %(inb,run,d)
rel='pc'
fil='0'

pfile=datafolder+"asc0_hmm_parameters_fil%s/p1_file" %(fil)
phalf=datafolder+"asc0_hmm_parameters_fil%s/phalf_file" %(fil)
datafile=datafolder+"asc0_mergedwin_fil%s/pw_%s.csv" %(fil,pair)
Afile="/mnt/diversity/divyaratan_popli/100arc/inbreeding/simulation_markov_inbreeding3/new_transitions1e6/transition_matrix_%s.csv" %(rel)
pairf1="/mnt/diversity/divyaratan_popli/100arc/inbreeding/simulation_markov_inbreeding3/hbd_hmm/inbred%s/run%s/downsample%s/filtered%s/asc0_res/gamma/pw_%s.csv" %(inb,run,d,fil,p1)
pairf2="/mnt/diversity/divyaratan_popli/100arc/inbreeding/simulation_markov_inbreeding3/hbd_hmm/inbred%s/run%s/downsample%s/filtered%s/asc0_res/gamma/pw_%s.csv" %(inb,run,d,fil,p2)
p1thresh=datafolder+"asc0_identicalmergedwin_fil1/pw_%s.csv" %(p1)
p2thresh=datafolder+"asc0_identicalmergedwin_fil1/pw_%s.csv" %(p2)
thresh=2
instates=[0,1,2]
inbr_states=3
inerror=0.001
upA=rel
"""


def getRelatable(filist, relatable, pairs, rels):

    cols=["pair", "relatedness","second_guess", "loglik_ratio", "withinDeg_second_guess", "withinDeg_ll"]
    df=pd.DataFrame(columns=cols)

    for i,fi in enumerate(filist):
        pair=pairs[i]
        lik=np.loadtxt(fi,dtype='float', delimiter = ",")
        data=pd.DataFrame(lik.reshape(-1,len(rels)), columns=rels)
        print(data)

        reldeg =	{
        "identical": "id",
        "siblings": "sib",
        "parent-child": "pc",
        "second": ["gr","hsib","avu"],
        "deg3": "deg3",
        "un": ["deg4","deg5", "un"],
        }


        b=data.sort_values(by=0, ascending=False, axis=1)

        df.loc[i,"relatedness"]=b.columns[0]
        df.loc[i,"pair"]=pair

        for r in reldeg:
            if b.columns[0] in reldeg[r]:
                temp=reldeg[r]

        if isinstance(temp, str):
            temp=[temp]

        b_otherdeg=b[list(set(b.columns) - set(temp))].sort_values(by=0, ascending=False, axis=1)
        df.loc[i,"second_guess"]=b_otherdeg.columns[0]
        df.loc[i,"loglik_ratio"]=b.iloc[0,0] - b_otherdeg.iloc[0,0]

        b_indeg=b[temp].sort_values(by=0, ascending=False, axis=1)
        if len(b_indeg.iloc[0]) > 1:
            df.loc[i,"withinDeg_second_guess"]=b_indeg.columns[1]
            df.loc[i,"withinDeg_ll"]=b_indeg.iloc[0,0] - b_indeg.iloc[0,1]
        else:
            df.loc[i,"withinDeg_second_guess"]=np.nan
            df.loc[i,"withinDeg_ll"]=np.nan

    with pd.option_context('display.max_rows', len(df.index), 'display.max_columns', len(df.columns)):
            df.to_csv(relatable, sep=',')



def IBDstates(gammafs, relatable, IBDfile,fil,dfolder):
    rel=pd.read_csv(relatable, sep=",", header=0,index_col=0)[['pair','relatedness']]

    pair_all=[]
    rel_all=[]
    k0prop=[]
    k0abs=[]
    k1prop=[]
    k1abs=[]
    k2prop=[]
    k2abs=[]
    IBDlen_all=[]
    IBDnum_all=[]
    for p in range(len(rel)):
        px=rel.loc[p,'pair']
        rx=rel.loc[p,'relatedness']
        gname="./filtered%s/relation_%s_file/res/gamma/pw_%s.csv" %(fil,rx,px)
        gamma=np.loadtxt(gname, dtype='float', delimiter = ",")
        gamma_sum=np.sum(gamma,0)/np.sum(gamma)

        rname="./filtered%s/relation_%s_file/res/pw_%s.csv" %(fil,rx,px)
        res=np.loadtxt(rname, dtype='float', delimiter = ",")
        res_count=np.array([np.sum(res==3/4), np.sum(res==1), np.sum(res==1/2)])/ len(res)

        IBD1=np.where(res==0.75)[0]
        IBDinfo=np.split(IBD1, np.where(np.diff(IBD1) != 1)[0]+1)
        IBDlen=len(IBD1)
        IBDnum=len(IBDinfo)
        pair_all.append(px)
        rel_all.append(rx)
        k0prop.append(gamma_sum[1])
        k1prop.append(gamma_sum[0])
        k2prop.append(gamma_sum[2])
        k0abs.append(res_count[1])
        k1abs.append(res_count[0])
        k2abs.append(res_count[2])
        IBDlen_all.append(IBDlen)
        IBDnum_all.append(IBDnum)

    df=pd.DataFrame({
        'pair': pair_all,
        'rel': rel_all,
        'k0_prop':k0prop,
        'k1_prop':k1prop,
        'k2_prop':k2prop,
        'k0_abs':k0abs,
        'k1_abs':k1abs,
        'k2_abs':k2abs,
        'IBD_len':IBDlen_all,
        'IBD_num':IBDnum_all
        })

    with pd.option_context('display.max_rows', len(df.index), 'display.max_columns', len(df.columns)):
        df.to_csv(IBDfile, sep=',')



def mergeRelatable(relf, mergef):


    rel=pd.read_csv(relf, sep=",", header=0,index_col=0)
    rel.loc[rel['relatedness']=='gr', 'relatedness'] = 'sec'
    rel.loc[rel['relatedness']=='hsib', 'relatedness'] = 'sec'
    rel.loc[rel['relatedness']=='avu', 'relatedness'] = 'sec'

    rel.loc[rel['relatedness']=='deg4', 'relatedness'] = 'un'
    rel.loc[rel['relatedness']=='deg5', 'relatedness'] = 'un'

    rel.loc[rel['relatedness']=='sec', 'withinDeg_ll'] = ''
    rel.loc[rel['relatedness']=='sec', 'withinDeg_second_guess'] = ''

    rel.loc[rel['relatedness']=='un', 'withinDeg_ll'] = ''
    rel.loc[rel['relatedness']=='un', 'withinDeg_second_guess'] = ''

    rel.loc[rel['second_guess']=='gr', 'second_guess'] = 'sec'
    rel.loc[rel['second_guess']=='hsib', 'second_guess'] = 'sec'
    rel.loc[rel['second_guess']=='avu', 'second_guess'] = 'sec'

    rel.loc[rel['second_guess']=='deg4', 'second_guess'] = 'un'
    rel.loc[rel['second_guess']=='deg5', 'second_guess'] = 'un'

    with pd.option_context('display.max_rows', len(rel.index), 'display.max_columns', len(rel.columns)):
        rel.to_csv(mergef, sep='\t')


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
