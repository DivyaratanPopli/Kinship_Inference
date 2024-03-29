
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
from scipy.stats import beta



def Betal(x,y):
    return gammaln(x)+gammaln(y) - gammaln(x+y)

def fact(x):
    return(gammaln(x+1))


def scaleB(B):
    scale=B.max(axis=0)
    Bscaled=B/scale
    return Bscaled, np.sum(np.log(scale))
    #return B, np.zeros(np.shape(B)[1])

def forward(A,B,pi,data):
    N= np.shape(A)[0]
    T= np.shape(data)[0]
    alpha=np.zeros((T,N))
    alpha_p=np.zeros((T,N))
    scale=np.zeros(T)

    alpha[0,:]= (pi) * B[:,0]
    alpha_p[0,:]= alpha[0,:]/sum(alpha[0,:])
    scale[0]=sum(alpha[0,:])

    for t in range(1,T):
        alpha[t,:]= (alpha_p[t-1,:] @ A) * B[:,t]
        alpha_p[t,:]= alpha[t,:]/sum(alpha[t,:])
        scale[t]=sum(alpha[t,:])

    return alpha_p, scale



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



def objective_old(x,gamma,data):
    if len(x)==6:
        a0=x[0]
        b0=x[1]
        a1=x[2]
        b1=x[3]
        a2=x[4]
        b2=x[5]
        cost= np.sum((fact(data[:,1]) - fact(data[:,0]) - fact(data[:,1]-data[:,0]) + Betal(data[:,0]+a0, data[:,1]-data[:,0]+b0) - Betal(a0,b0))*gamma[0,:]) + \
              np.sum((fact(data[:,1]) - fact(data[:,0]) - fact(data[:,1]-data[:,0]) + Betal(data[:,0]+a1, data[:,1]-data[:,0]+b1) - Betal(a1,b1))*gamma[1,:]) + \
              np.sum((fact(data[:,1]) - fact(data[:,0]) - fact(data[:,1]-data[:,0]) + Betal(data[:,0]+a2, data[:,1]-data[:,0]+b2) - Betal(a2,b2))*gamma[2,:])

    elif len(x)==4:
        a0=x[0]
        b0=x[1]
        a1=x[2]
        b1=x[3]

        cost= np.sum((fact(data[:,1]) - fact(data[:,0]) - fact(data[:,1]-data[:,0]) + Betal(data[:,0]+a0, data[:,1]-data[:,0]+b0) - Betal(a0,b0))*gamma[0,:]) + \
              np.sum((fact(data[:,1]) - fact(data[:,0]) - fact(data[:,1]-data[:,0]) + Betal(data[:,0]+a1, data[:,1]-data[:,0]+b1) - Betal(a1,b1))*gamma[1,:])


    return (-1)*cost




def objective(a, meanp, k, data):
    cost=0

    b=a * (1-meanp)/meanp

    for w in range(len(data)):

        ll=(fact(data[w,1]) - fact(data[w,0]) - fact(data[w,1]-data[w,0]) + \
             Betal(data[w,0]+a, data[w,1]-data[w,0]+b) - Betal(a,b))


        Wll= ll * k[w]

        cost= cost + Wll

    return (-1)*cost





def makeB(data,up_a0,up_b0,up_a1,up_b1,M,T,hd):


    B = np.zeros((M,T))
    B[0,:]= np.exp(fact(data[:,1]) - fact(data[:,0]) - fact(data[:,1]-data[:,0]) + Betal(data[:,0]+up_a0, data[:,1]-data[:,0]+up_b0) - Betal(up_a0,up_b0))
    B[1,:]= np.exp(fact(data[:,1]) - fact(data[:,0]) - fact(data[:,1]-data[:,0]) + Betal(data[:,0]+up_a1, data[:,1]-data[:,0]+up_b1) - Betal(up_a1,up_b1))

    B[0,hd]=0
    B[1,hd]=1
    return B


def expectation(data, A, B, pos, chrm_no, log_bscale, pi):

    M = A.shape[0]
    T = len(data)
    alpha=np.zeros((T,M))
    beta=np.zeros((T,M))
    scale=np.zeros(T)




    for chrm in range(chrm_no):
        data_chr=data[pos[chrm]:pos[chrm+1],:]

        alpha[pos[chrm]:pos[chrm+1],:],scale[pos[chrm]:pos[chrm+1]] = forward(A,B[:,pos[chrm]:pos[chrm+1]],pi,data_chr)
        beta[pos[chrm]:pos[chrm+1],:] = backward(A,B[:,pos[chrm]:pos[chrm+1]],data_chr,scale[pos[chrm]:pos[chrm+1]])

    post= alpha*beta


    li=np.sum(np.log(scale)) + log_bscale


    assert isclose(sum(abs(np.sum(post,axis=1) - np.ones(np.shape(B)[1]))), 0, abs_tol=1e-7), "sum of alpha and beta is not 1"
    assert isclose(sum(abs(np.sum(A,axis=1) - np.ones(np.shape(A)[0]))), 0, abs_tol=1e-7), "sum of transition matrix columns is not 1"

    #pi_new=np.mean(post[pos[:-1],:],0)
    pi_new=matrix_power(A,10000)[0]
    return post.T, li, alpha, beta, scale, pi_new




def transition(alpha1, beta1, n1, gamma1, A, B1, pos1):
    pos=pos1[1:-1]
    B=np.vsplit(B1.T, pos)
    gamma=np.vsplit(gamma1.T, pos)
    alpha=np.vsplit(alpha1, pos)
    beta=np.vsplit(beta1, pos)
    n=np.hsplit(n1, pos)



    A_new = np.zeros_like(A)
    n_states = A.shape[0]

    for i in range(n_states):
        for j in range(n_states):
            for a, b, e, n_ in zip(alpha, beta, B, n):
                A_new[i, j] += np.sum(
                    a[:-1, i] * A[i, j] * b[1:, j] * e[1:, j] / n_[1:]
                )

    gamma_sum = np.sum([np.sum(g[:-1], 0) for g in gamma], 0)
    A_new /= gamma_sum[:, np.newaxis]
    assert np.allclose(np.sum(A_new, 1), 1)

    return A_new

def emission(gamma, data, p1avg, M, x0, highdiff):

    T = len(data)



    # optimize
    bd = (0.01,100000)
    mean_inb,mean_id=p1avg

    solution_inb = minimize_scalar(objective,x0[0],args=(mean_inb,gamma[0,:],data),method='Bounded',\
                        bounds=bd)
    solution_id = minimize_scalar(objective,x0[2],args=(mean_id,gamma[1,:],data),method='Bounded',\
                        bounds=bd)

    a_inb=solution_inb.x
    a_id=solution_id.x
    b_inb= a_inb * (1-mean_inb)/mean_inb
    b_id= a_id * (1-mean_id)/mean_id

    x=[a_inb,b_inb,a_id,b_id]

    up_p=np.zeros((M))


    up_a0=x[0]
    up_b0=x[1]
    up_a1=x[2]
    up_b1=x[3]


    up_p[0]=up_a0/(up_a0+up_b0)
    up_p[1]=up_a1/(up_a1+up_b1)



    B=makeB(data=data, up_a0=up_a0, up_b0=up_b0, up_a1=up_a1, up_b1=up_b1, M=M, T=T, hd=highdiff)
    [B, log_bscale]=scaleB(B)

    return B, up_p, log_bscale,x


def bnd_calc(mean_p, dist, propvar):

    r=(1-mean_p)/mean_p
    t=propvar*dist*dist

    low_bnd=-1 * (t + (t*r*r) + r*((2*t) - 1))/(t+ (3*t*r) + (3*r*r*t) + (r*r*r*t))
    return low_bnd

def baum_welch(data, A, B, pos, p1avg, log_bscale, x0, highdiff, max_iter=1000):

    chrm_no=len(pos)-1
    n_iter=0
    diff=1
    last_lik=0



    #print("initial beta parameters are:",x0)

    [mean_inb,mean_id]=[p1avg[0],p1avg[1]]

    di=mean_id-mean_inb


    bnd_inb= bnd_calc(mean_p=mean_inb, dist=di, propvar=1)
    bnd_id= bnd_calc(mean_p=mean_id, dist=di, propvar=1)
    bnds=[bnd_inb,bnd_id]
    #I have not used any bounds yet
    pi=[0.5,0.5]

    while diff > 0.000001 and n_iter<max_iter:



        [gamma, li, alpha, beta, scale, pi]=expectation(data=data, A=A, B=B, pos=pos, chrm_no=chrm_no, log_bscale=log_bscale, pi=pi)
        #print("alpha",alpha[:20,:],n_iter)
        A=transition(alpha1=alpha, beta1=beta, n1=scale, gamma1=gamma, A=A, B1=B, pos1=pos)


        [B, up_p, log_bscale,x1]=emission(gamma=gamma, data=data, p1avg=p1avg, M=A.shape[0], x0=x0, highdiff=highdiff)
        x0=x1

        if n_iter !=0:
            diff= li - last_lik
            last_lik=li
        elif n_iter==0:
            diff= last_lik - li
            last_lik=li
        n_iter=n_iter+1
        #print("diff=",diff)
        #print("up_p=", up_p)
        #print("x=",x0)
        #print("lik=",li)
        #print("A=",A)


    return gamma, A , B, li, x1



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
    #S[S==2]=1/2
    #S[S==0]=3/4
    return S



def hmm(dfile, tfile, ind, chrmf, p1file, upA, targets, mult):

    ind1=np.where(np.array(targets)==ind)[0][0]
    diff1= np.loadtxt(dfile, dtype='float', delimiter = ",")[:,ind1]
    total1=np.loadtxt(tfile, dtype='float', delimiter = ",")[:,ind1]

    chrm1=np.loadtxt(chrmf, dtype='float', delimiter = ",")

    totalwin=len(total1)
    gudwin=np.where(total1>0)

    diff=diff1[gudwin]
    total=total1[gudwin]
    chrm=chrm1[gudwin]

    win=len(total)


    with open(p1file,"r") as f:
        p_1=float(f.read())

    p_12=p_1/2
    p_in=p_12/2
    initial_p=np.array([p_in, p_12])

    highdiff=np.where(diff/total>p_1*mult)


    data=np.stack([diff,total]).T


    pos=np.where(chrm[:-1] != chrm[1:])[0]+1
    pos=np.append(pos,np.shape(total)[0])
    pos=np.insert(pos, 0, 0, axis=0)


    A=np.array([[0.8,0.2],[0.2,0.8]])





    Bu = np.zeros((np.shape(A)[0],win)) #Emission probability
    [b0, b1]=[1000,1000]
    [a0,a1]=[b0*initial_p[0]/(1-initial_p[0]), b1*initial_p[1]/(1-initial_p[1])]
    x0 = [a0,b0,a1,b1]

    Bu=makeB(data=data, up_a0=x0[0], up_b0=x0[1], up_a1=x0[2], up_b1=x0[3], M=np.shape(A)[0], T=len(data), hd=highdiff)


    [B, log_bscale]= scaleB(Bu)

    np.argwhere(np.isnan(B))

    #print("initial A=",A)

    if win>40:              #if number of windows with nonzero count
        gamma,A,B,lik,betas= baum_welch(data=data, A=A, B=B, pos=pos, p1avg=initial_p, log_bscale=log_bscale, x0=x0, highdiff=highdiff)
    else:
        gamma=np.ones([2,win]) * 9
        lik='9'
    #res=viterbi(data, A, B, pi)

    #resall=np.ones(totalwin)*9
    #resall[gudwin]=res
    #np.savetxt(fname=resfile, X=resall,delimiter=',')
    #observations1['count']=1
    #observations1['dis']=resall
    gammall=np.ones((totalwin,2))*9
    gammall[gudwin,:]=gamma.T

    dfg=pd.DataFrame({'chrom':chrm1,
                        'count': total1,
                        'g_noin': gammall[:,1]})

    harr=np.asarray(dfg['g_noin'])
    goodprop=np.sum((harr>=0.95) & (harr<1.1))/np.sum(harr<=1.1)
    if goodprop<0.5:
        dfg['g_noin']=9



    return gamma,A,B,lik,betas,highdiff


"""
#p1="Chagyrskaya07"

cov=0.2
ind=2
m=1

datafolder="/mnt/diversity/divyaratan_popli/100arc/inbreeding/fastsim_gz/"

difffile=datafolder+"contam0/inbred1/run57/coverage%s/asc0/inputMode_idhapProbs_fil0/id_remHighdiv_diff.csv.gz" %(cov)
totalfile=datafolder+"contam0/inbred1/run57/coverage%s/asc0/inputMode_idhapProbs_fil0/id_remHighdiv_total.csv.gz" %(cov)

p1file=datafolder+"contam0/inbred1/run57/coverage%s/asc0/inputMode_hapProbs_fil0/hmm_parameters/p1_file.gz" %(cov)
chrmf= datafolder + "contam0/inbred1/run57/coverage%s/asc0/inputMode_idhapProbs_fil0/id_chrm.csv.gz" %(cov)
targets=list(range(17))
gamma,A,B,lik,betas,hd1=hmm(dfile=difffile, tfile=totalfile, ind=ind, chrmf=chrmf,
                            p1file=p1file, upA="un", targets=targets, mult=m)


"""

#betasnf=[1.233898231777093, 85.65018006593586, 5.807562635082375, 198.66016157044342]
#betasf=[1.1538670587410393, 80.09487233885226, 5.842327638089113, 199.84937321537603]



def plotbetas(dfile, tfile, ind, chrmf, p1file, upA, targets, outplot):
    ind=int(ind)
    mall=[100,1]
    x = np.linspace(0,0.1,1000)
    fig, ax = plt.subplots(nrows=1, ncols=2)
    for i in [0,1]:
        m=mall[i]

        gamma,A,B,lik,betas,hd1=hmm(dfile=dfile, tfile=tfile, ind=ind, chrmf=chrmf,
                                    p1file=p1file, upA="un", targets=targets, mult=m)

        ind1=np.where(np.array(targets)==ind)[0][0]
        diff= np.loadtxt(dfile, dtype='float', delimiter = ",")[:,ind1]
        total=np.loadtxt(tfile, dtype='float', delimiter = ",")[:,ind1]



        a=betas[0]
        b=betas[1]
        ax[i].plot(x, beta.pdf(x, a, b)/500,'b-', lw=3, alpha=0.6, label="$Y_w$=4")

        a=betas[2]
        b=betas[3]
        ax[i].plot(x, beta.pdf(x, a, b)/500,'r-', lw=3, alpha=0.6, label="$Y_w$=2")
        ax[i].hist(diff/total,alpha=0.3, weights=np.ones(len(diff)) / len(diff))

        ax[i].set_ylim(0,1)
        ax[i].set_xlim(0,0.04)

    ax[0].set_title('Without Constraint')
    ax[1].set_title('With Constraint')

    ax[1].legend(loc="upper right")
    fig.text(0.5, 0.04, 'Proportion of differences', ha='center')
    fig.text(0.04, 0.5, 'Frequency', va='center', rotation='vertical')
    plt.savefig(outplot, bbox_inches='tight')
