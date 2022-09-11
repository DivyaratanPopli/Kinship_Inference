
import scipy
import numpy as np
import random
import pandas as pd
from scipy.stats import binom
from math import isclose
from scipy.special import gammaln
from scipy import optimize
import math
from scipy.special import loggamma
from scipy.special import beta as Betaf
from scipy.optimize import minimize_scalar
import pylab
from numpy.linalg import matrix_power
from scipy.special import logsumexp
import numba
from scipy.stats import beta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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

def baum_welch(data, hbd, A, B, pos, p1avg, update_transition, inbr, x0, mode, max_iter=1000):
    #print(data)
    chrm_no=len(pos)-1
    n_iter=0
    diff=1
    last_lik=0


    #print("initial beta parameters are:",x0)

    [mean_pc,mean_un,mean_id,mean_inb]=[p1avg[0,0],p1avg[0,1],p1avg[0,2],p1avg[2,2]]

    di=(mean_un-mean_pc)

    if mode=='bnds':
        bnd_pc= bnd_calc(mean_p=mean_pc, dist=di, propvar=3/4)
        bnd_un= bnd_calc(mean_p=mean_un, dist=di, propvar=1)
        bnd_id= bnd_calc(mean_p=mean_id, dist=di, propvar=1/2)
        bnd_inb= bnd_calc(mean_p=mean_inb, dist=di, propvar=1)
    else:
        bnd_pc= float(mode)
        bnd_un= float(mode)
        bnd_id= float(mode)
        bnd_inb= float(mode)

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
    return gamma, A , B2final, up_p, li, pi,x0


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

def hmm(difffile, targets, totalfile, listind, listf, pfile, Afile, upA, pairf1, pairf2, instates, p1thresh, thresh, mode):


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



        pi= np.array([1/3,1/3,1/3])

        #A=A1
        A=np.array(pd.read_csv(Afile,sep=',',header=None,index_col=False))

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



        gamma,A,B,up_p,lik, pi,xnew= baum_welch(data=data, hbd=hbd, A=A, B=B, pos=pos, p1avg=initial_p, update_transition=upA, inbr=inbr_states, x0=x0, mode=mode)
        res=viterbi(data, A, B, pi)

        B_all2=np.ones(len(total))*-9
        B_all2[gudwin]=B[2,:]

    else:
        res=np.ones(len(data))*-9
        lik=float('-inf')
        gamma=np.ones([len(data),3])*-9

    #np.savetxt(fname=resfile, X=res,delimiter=',')
    #np.savetxt(fname=gammafile, X=gamma.T,delimiter=',')

    #with open(likfile,'w') as f:
    #    f.write(str(lik))

    #print(lik)
    return res,lik,gamma,B_all2,xnew



"""

cov=0.1
cnt="0"
inb="0"
p1="9"
p2="10"
pair=p1+"_"+p2
run=1
datafolder="/mnt/diversity/divyaratan_popli/100arc/inbreeding/fastsim_gz/"
rel='pc'
fil='0'

pfile=datafolder+"contam%s/inbred%s/run%s/coverage%s/asc0/inputMode_hapProbs_fil0/hmm_parameters/p1_file.gz" %(cnt,inb,run,cov)
dfile=datafolder+"contam%s/inbred%s/run%s/coverage%s/asc0/inputMode_hapProbs_fil0/merged_wind.csv.gz" %(cnt,inb,run,cov)
tfile=datafolder+"contam%s/inbred%s/run%s/coverage%s/asc0/inputMode_hapProbs_fil0/merged_wint.csv.gz" %(cnt,inb,run,cov)

Afile="/mnt/diversity/divyaratan_popli/100arc/inbreeding/fastsim_gz/transitions/transition_matrix_%s.csv" %(rel)
pairf1=datafolder+"hbd_hmm/contam%s/inbred%s/run%s/coverage%s/filtered0/asc0_res/gamma_hapProbs/pw_%s.csv.gz" %(cnt,inb,run,cov,p1)
pairf2=datafolder+"hbd_hmm/contam%s/inbred%s/run%s/coverage%s/filtered0/asc0_res/gamma_hapProbs/pw_%s.csv.gz" %(cnt,inb,run,cov,p2)
p1thresh=datafolder+"contam%s/inbred%s/run%s/coverage%s/asc0/inputMode_idhapProbs_fil1/id_wint.csv.gz" %(cnt,inb,run,cov)
thresh=10
instates=[0,1,2]
inbr_states=3
upA=rel
libraries=list(range(17))
listf=[]
for i, l1 in enumerate(libraries):
    for l2 in libraries[(i+1):]:
        s = '%s_%s' %(str(l1).strip("['']"), str(l2).strip("['']"))
        listf.append(s)




res1,lik1,gamma1,x1=hmm(difffile=dfile, targets=np.array(libraries), totalfile=tfile, listind=pair, listf=np.array(listf), pfile=pfile, Afile=Afile, resfile='resfile', likfile='likfile', gammafile='gammafile', upA=upA, pairf1=pairf1, pairf2=pairf2, instates=instates, p1thresh=p1thresh, thresh=thresh)

#looping

lik11=[]
lik12=[]

for r in range(1,61):
    for rel in ['id','sib']:
        cov=0.2
        cnt="0"
        inb="1"
        p1="0"
        p2="15"
        pair=p1+"_"+p2
        run=r
        datafolder="/mnt/diversity/divyaratan_popli/100arc/inbreeding/fastsim_gz/"
        #rel='sib'
        fil='0'

        pfile=datafolder+"contam%s/inbred%s/run%s/coverage%s/asc0/inputMode_hapProbs_fil0/hmm_parameters/p1_file.gz" %(cnt,inb,run,cov)
        dfile=datafolder+"contam%s/inbred%s/run%s/coverage%s/asc0/inputMode_hapProbs_fil0/remHighdiv_diff.csv.gz" %(cnt,inb,run,cov)
        tfile=datafolder+"contam%s/inbred%s/run%s/coverage%s/asc0/inputMode_hapProbs_fil0/remHighdiv_total.csv.gz" %(cnt,inb,run,cov)

        Afile="/mnt/diversity/divyaratan_popli/100arc/inbreeding/fastsim_gz/transitions/transition_matrix_%s.csv" %(rel)
        pairf1=datafolder+"hbd_hmm/contam%s/inbred%s/run%s/coverage%s/filtered0/asc0_res/gamma_hapProbs/pw_%s.csv.gz" %(cnt,inb,run,cov,p1)
        pairf2=datafolder+"hbd_hmm/contam%s/inbred%s/run%s/coverage%s/filtered0/asc0_res/gamma_hapProbs/pw_%s.csv.gz" %(cnt,inb,run,cov,p2)
        p1thresh=datafolder+"contam%s/inbred%s/run%s/coverage%s/asc0/inputMode_idhapProbs_fil1/id_remHighdiv_total.csv.gz" %(cnt,inb,run,cov)
        thresh=10
        instates=[0,1,2]
        inbr_states=3
        upA=rel
        libraries=list(range(17))
        listf=[]
        for i, l1 in enumerate(libraries):
            for l2 in libraries[(i+1):]:
                s = '%s_%s' %(str(l1).strip("['']"), str(l2).strip("['']"))
                listf.append(s)




        res1,lik1,gamma1,B,xnew=hmm(difffile=dfile, targets=np.array(libraries), totalfile=tfile, listind=pair, listf=np.array(listf), pfile=pfile, Afile=Afile, resfile='resfile', likfile='likfile', gammafile='gammafile', upA=upA, pairf1=pairf1, pairf2=pairf2, instates=instates, p1thresh=p1thresh, thresh=thresh)
        if rel=='id':
            lik11.append(lik1)
        elif rel=='sib':
            lik12.append(lik1)




ind=np.where(np.array(listf)==pair)[0][0]



diff= np.loadtxt(dfile,dtype='float', delimiter = ",")[:,ind]
total= np.loadtxt(tfile,dtype='float', delimiter = ",")[:,ind]

plt.plot(diff/total)
plt.plot(B[2,:]/40)




betas=xnew
x = np.linspace(0,0.06,1000)

a=betas[0]
b=betas[1]
r1 = beta.rvs(a, b, size=1000)
plt.plot(x, beta.pdf(x, a, b)/1000,'b-', lw=6, alpha=0.6, label='beta pdf')

a=betas[2]
b=betas[3]
r1 = beta.rvs(a, b, size=1000)
plt.plot(x, beta.pdf(x, a, b)/1000,'r-', lw=6, alpha=0.6, label='beta pdf')
#plt.hist(diff/total,alpha=0.3)


a=betas[4]
b=betas[5]
r1 = beta.rvs(a, b, size=1000)
plt.plot(x, beta.pdf(x, a, b)/1000,'g-', lw=6, alpha=0.6, label='beta pdf')
#plt.hist(diff/total,alpha=0.3)

a=betas[6]
b=betas[7]
r1 = beta.rvs(a, b, size=1000)
plt.plot(x, beta.pdf(x, a, b)/1000,'y-', lw=6, alpha=0.6, label='beta pdf')
plt.hist(diff/total,alpha=0.3,weights=np.ones(len(diff)) / len(diff))
plt.ylim(0,1)

"""

def plot_betas(dfile, tfile, pair, thresh, rel, libraries, listf, Afile, pfile, pairf1, pairf2, p1thresh, instates, B1file, B2file, xnew1file, xnew2file, res1file, res2file, betaplot):
    #cov=0.2
    #cnt="0"
    #inb="1"
    #p1="0"
    #p2="15"
    #pair=p1+"_"+p2
    #run=r
    #datafolder="/mnt/diversity/divyaratan_popli/100arc/inbreeding/fastsim_gz/"
    #rel='sib'

    #pfile=datafolder+"contam%s/inbred%s/run%s/coverage%s/asc0/inputMode_hapProbs_fil0/hmm_parameters/p1_file.gz" %(cnt,inb,run,cov)
    #dfile=datafolder+"contam%s/inbred%s/run%s/coverage%s/asc0/inputMode_hapProbs_fil0/remHighdiv_diff.csv.gz" %(cnt,inb,run,cov)
    #tfile=datafolder+"contam%s/inbred%s/run%s/coverage%s/asc0/inputMode_hapProbs_fil0/remHighdiv_total.csv.gz" %(cnt,inb,run,cov)

    #Afile="/mnt/diversity/divyaratan_popli/100arc/inbreeding/fastsim_gz/transitions/transition_matrix_%s.csv" %(rel)
    #pairf1=datafolder+"hbd_hmm/contam%s/inbred%s/run%s/coverage%s/filtered0/asc0_res/gamma_hapProbs/pw_%s.csv.gz" %(cnt,inb,run,cov,p1)
    #pairf2=datafolder+"hbd_hmm/contam%s/inbred%s/run%s/coverage%s/filtered0/asc0_res/gamma_hapProbs/pw_%s.csv.gz" %(cnt,inb,run,cov,p2)
    #p1thresh=datafolder+"contam%s/inbred%s/run%s/coverage%s/asc0/inputMode_idhapProbs_fil1/id_remHighdiv_total.csv.gz" %(cnt,inb,run,cov)
    #thresh=10
    #instates=[0,1,2]
    #upA=rel
    #libraries=list(range(17))
    mode='0.1'
    res1,lik1,gamma1,B1,xnew1=hmm(difffile=dfile, targets=np.array(libraries), totalfile=tfile, listind=pair, listf=np.array(listf), pfile=pfile, Afile=Afile, upA=rel, pairf1=pairf1, pairf2=pairf2, instates=instates, p1thresh=p1thresh, thresh=thresh, mode=mode)
    mode='bnds'
    res2,lik2,gamma2,B2,xnew2=hmm(difffile=dfile, targets=np.array(libraries), totalfile=tfile, listind=pair, listf=np.array(listf), pfile=pfile, Afile=Afile, upA=rel, pairf1=pairf1, pairf2=pairf2, instates=instates, p1thresh=p1thresh, thresh=thresh, mode=mode)

    xall=[xnew1,xnew2]
    fig, ax = plt.subplots(nrows=1, ncols=2)
    ind=np.where(np.array(listf)==pair)[0][0]

    diff= np.loadtxt(dfile,dtype='float', delimiter = ",")[:,ind]
    total= np.loadtxt(tfile,dtype='float', delimiter = ",")[:,ind]
    for i in [0,1]:
        betas=xall[i]
        x = np.linspace(0,0.15,1000)

        #no. of points in each bin=1000*0.08/0.15
        #I normalize with this value to make it comparable to the prop. of diff. histogram which is also normalized by the
        #no. of total windows

        a=betas[2]
        b=betas[3]
        ax[i].plot(x, beta.pdf(x, a, b)/533.33,'r-', lw=3, alpha=0.6, label="i=0")
        #plt.hist(diff/total,alpha=0.3)

        a=betas[0]
        b=betas[1]
        ax[i].plot(x, beta.pdf(x, a, b)/533.33,'b-', lw=3, alpha=0.6, label="i=1")

        a=betas[4]
        b=betas[5]
        ax[i].plot(x, beta.pdf(x, a, b)/533.33,'g-', lw=3, alpha=0.6, label="i=2")
        #plt.hist(diff/total,alpha=0.3)

        a=betas[6]
        b=betas[7]

        ax[i].plot(x, beta.pdf(x, a, b)/533.33,'y-', lw=3, alpha=0.6, label="i=4")
        ax[i].hist(diff/total,alpha=0.3,weights=np.ones(len(total)) / len(total))
        ax[i].set_ylim(0,0.4)

    ax[0].set_title('Without Constraint')
    ax[1].set_title('With Constraint')

    ax[1].legend(loc="upper right")
    fig.text(0.5, 0.04, 'Frequency', ha='center')
    fig.text(0.04, 0.5, 'Proportion of differences', va='center', rotation='vertical')
    plt.savefig(betaplot, bbox_inches='tight')

    np.savetxt(fname=B1file, X=B1,delimiter=',')
    np.savetxt(fname=B2file, X=B2,delimiter=',')
    np.savetxt(fname=xnew1file, X=xnew1,delimiter=',')
    np.savetxt(fname=xnew2file, X=xnew2,delimiter=',')
    np.savetxt(fname=res1file, X=res1,delimiter=',')
    np.savetxt(fname=res2file, X=res2,delimiter=',')
