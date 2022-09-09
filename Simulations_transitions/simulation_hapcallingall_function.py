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

import copy
from ast import literal_eval
import shutil
import scipy.special as scp
##########################################################start with simulation for relatedness distn


def get_recomb(recombination_rate_child,length):

        lam1=recombination_rate_child*length
        recomb_no1=np.random.poisson(lam1)
        recomb_pos1=np.zeros(recomb_no1)
        recomb_pos1=np.random.uniform(0,length,recomb_no1)
        recomb_pos1=np.sort(recomb_pos1)

        recomb_interval1=np.zeros(recomb_no1+2)
        recomb_interval1[0]=0
        recomb_interval1[-1]=length
        recomb_interval1[1:recomb_no1+1]=recomb_pos1

        lam2=recombination_rate_child*length
        recomb_no2=np.random.poisson(lam2)
        recomb_pos2=np.zeros(recomb_no2)
        recomb_pos2=np.random.uniform(0,length,recomb_no2)
        recomb_pos2=np.sort(recomb_pos2)


        recomb_interval2=np.zeros(recomb_no2+2)
        recomb_interval2[0]=0
        recomb_interval2[-1]=length
        recomb_interval2[1:recomb_no2+1]=recomb_pos2

        return recomb_no1, recomb_interval1, recomb_no2, recomb_interval2




def make_baby(gen, pos, parent1, parent2, recomb_no1, recomb_interval1, recomb_no2, recomb_interval2, x1, x2):
    #recombination_rate_child=1e-8     #increasing recombinations for artificial kids 2e-6, 4e-7 for len=1e7


    c1=np.zeros(np.shape(gen)[0])
    c1alt=np.zeros(np.shape(gen)[0])
    for i in range(recomb_no1+1):
        xx=np.where(np.logical_and(pos>=recomb_interval1[i], pos<recomb_interval1[i+1]))
        if i%2==0:
            c1[xx] = gen[xx,parent1[0]]
            c1alt[xx]=gen[xx,parent1[1]]
        elif i%2==1:
            c1[xx] = gen[xx,parent1[1]]
            c1alt[xx]=gen[xx,parent1[0]]

    c2=np.zeros(np.shape(gen)[0])
    c2alt=np.zeros(np.shape(gen)[0])
    for i in range(recomb_no2+1):
        xx=np.where(np.logical_and(pos>=recomb_interval2[i], pos<=recomb_interval2[i+1]))
        if i%2==0:
            c2[xx] = gen[xx,parent2[0]]
            c2alt[xx]=gen[xx,parent2[1]]
        elif i%2==1:
            c2[xx] = gen[xx,parent2[1]]
            c2alt[xx] = gen[xx,parent2[0]]

    chrmPair1=[c1,c1alt]
    chrmPair2=[c2,c2alt]


    return chrmPair1[x1], chrmPair2[x2], x1, x2, recomb_interval1, recomb_interval2


def make_inbred(len_file,a00,a11):

    l=len_file

    a=a00  #hbd
    b=a11  #non hbd

    trans=np.array([[a,1-a],[1-b,b]])
    pi=[1/2,1/2]
    states=[0,1]

    seq=np.zeros(l)

    seq[0]= np.random.binomial(1, 0.5)

    for i in range(1,l):
        p=trans[int(seq[i-1]),1]
        seq[i]=np.random.binomial(1, p)

    return seq


def simulate_seq(unrelated,
                 related,
                 identical,
                 others,
                 out_gen,
                 pos_file,
                 chrm,
                 RID,
                 Ne=1e3,
                 length=1e8,
                 recombination_rate=1e-8,
                 mutation_rate=1e-8,
                 recombination_rate_child=1e-8
                 ):

    chag=0
    altai=1
    vindija=2
    denisova=3
    human=4

    population_configurations = [
            msp.PopulationConfiguration(
                initial_size=Ne, growth_rate=0),
            msp.PopulationConfiguration(
                initial_size=Ne, growth_rate=0),
            msp.PopulationConfiguration(
                initial_size=Ne, growth_rate=0),
            msp.PopulationConfiguration(
                initial_size=Ne, growth_rate=0),
            msp.PopulationConfiguration(
                initial_size=10*Ne, growth_rate=0)]

    demo = [
            msp.MassMigration(time=3500,source=2,destination=0,proportion=1.0),
            msp.MassMigration(time=4500,source=0,destination=1,proportion=1.0),
            msp.MassMigration(time=14000,source=1,destination=3,proportion=1.0),
            msp.MassMigration(time=20000,source=3,destination=4,proportion=1.0)]


    samples = [msp.Sample(population=0, time=2500) ]*unrelated + [msp.Sample(population=1, time=4000) ]*int(others/8) + [msp.Sample(population=2, time=2000) ]*int(others/8) + [msp.Sample(population=3, time=2500) ]*int(others/8) + [msp.Sample(population=4, time=0) ]*int(others-(others*3/8))

    tree_sequence = msp.simulate(
            length=length, recombination_rate=recombination_rate,
            mutation_rate=mutation_rate,population_configurations=population_configurations,
            demographic_events=demo,samples=samples)

    for j, variant in enumerate(tree_sequence.variants()):

        x=variant.site.id, variant.site.position,
        variant.alleles, variant.genotypes

    no_sites=variant.site.id+1
    gen= np.zeros((no_sites, (unrelated+related+identical+others)))
    pos= np.zeros(no_sites)

    for i,variant in enumerate(tree_sequence.variants()):
        gen[i,0:(unrelated)]= variant.genotypes[0:(unrelated)]
        gen[i,(unrelated+related+identical):(unrelated+related+identical+others)]= variant.genotypes[unrelated:(unrelated+others)]
        pos[i]= variant.site.position

    Dip=pd.DataFrame(
                data=gen,
                index=np.array(range(0,np.shape(gen)[0])),
                columns=np.array(range(0,np.shape(gen)[1]))
                )

    with pd.option_context('display.max_rows', len(Dip.index), 'display.max_columns', len(Dip.columns)):
                 Dip.to_csv(out_gen, sep=',')

    np.savetxt(fname=pos_file, X=pos, delimiter=",")



def postSim(genfile,
                 posfile,
                 unrelated,
                 related,
                 identical,
                 others,
                 diploid_sim,
                 out_hbdall,
                 chrm,
                 RID,
                 Ne=1e3,
                 length=1e8,
                 recombination_rate=1e-8,
                 mutation_rate=1e-8,
                 recombination_rate_child=1e-8):

    chrm_no='chrm' + chrm + '/'

    df_gen=pd.read_csv(genfile, sep=",", header=0,index_col=0)
    gen=np.asarray(df_gen)
    gen_i=copy.copy(gen)

    pos=np.loadtxt(posfile,dtype='float', delimiter = ",")

    hbdall=np.ones([len(gen),int((unrelated+related+identical)/2)])
    #hbdall_i=np.ones([len(gen),int((unrelated+related+identical)/2)])

    #introducing inbreeding
    gen_un=copy.copy(gen_i[:,0:(unrelated)])


    l=len(gen_un)


    #add hbd state to identical individuals
    #hbdall_i[:,int((unrelated+related)/2):int((unrelated+related+identical)/2)]=hbdall_i[:,:int((identical)/2)]

    col=list(range(int((unrelated+related+identical)/2)))


    #gen contains the chromosomes of all unrelated people, gen_i contains chromosomes with inbreeding

    #now making artificial chromosomes of four children by mating, also 1 grand kid, and 1 (3rd) relative

    p1=0,1    #parent1
    p2=2,3    #parent2
    p3=4,5    #parent3
    p4=6,7    #parent4
    p5=8,9
    p6=10,11
    p7=12,13

    c1_12=unrelated,unrelated+1    #child1 from mating of p1 and p2
    c2_23=unrelated+2,unrelated+3
    c3_23=unrelated+4,unrelated+5
    c4_45=unrelated+6,unrelated+7
    gc1_1011=unrelated+8,unrelated+9
    ggc1_512=unrelated+10,unrelated+11
    gggc1_613=unrelated+12,unrelated+13



    #getting recomb points
    [no1_c1_12, interv1_c1_12, no2_c1_12, interv2_c1_12]=get_recomb(length=length, recombination_rate_child=recombination_rate_child)
    [no1_c2_23, interv1_c2_23, no2_c2_23, interv2_c2_23]=get_recomb(length=length, recombination_rate_child=recombination_rate_child)
    [no1_c3_23, interv1_c3_23, no2_c3_23, interv2_c3_23]=get_recomb(length=length, recombination_rate_child=recombination_rate_child)
    [no1_c4_45, interv1_c4_45, no2_c4_45, interv2_c4_45]=get_recomb(length=length, recombination_rate_child=recombination_rate_child)
    [no1_gc1_1011, interv1_gc1_1011, no2_gc1_1011, interv2_gc1_1011]=get_recomb(length=length, recombination_rate_child=recombination_rate_child)
    [no1_ggc1_512, interv1_ggc1_512, no2_ggc1_512, interv2_ggc1_512]=get_recomb(length=length, recombination_rate_child=recombination_rate_child)
    [no1_gggc1_613, interv1_gggc1_613, no2_gggc1_613, interv2_gggc1_613]=get_recomb(length=length, recombination_rate_child=recombination_rate_child)

    #deciding randomly which parent will give the first fragment

    #x1=np.random.binomial(1,0.5)
    #x2=np.random.binomial(1,0.5)

    #rec1_12_1 means recombination even in child 1's chromosome from p1 and p2 in chromosome1
    gen[:,c1_12[0]], gen[:,c1_12[1]], x1_12_1,x1_12_2, rec1_12_1, rec1_12_2= make_baby(gen=gen, pos=pos, parent1=p1, parent2=p2,recomb_no1=no1_c1_12, recomb_interval1=interv1_c1_12, recomb_no2=no2_c1_12, recomb_interval2=interv2_c1_12, x1=np.random.binomial(1,0.5), x2=np.random.binomial(1,0.5))
    gen[:,c2_23[0]], gen[:,c2_23[1]], x2_23_1,x2_23_2, rec2_23_1, rec2_23_2= make_baby(gen=gen, pos=pos, parent1=p2, parent2=p3,recomb_no1=no1_c2_23, recomb_interval1=interv1_c2_23, recomb_no2=no2_c2_23, recomb_interval2=interv2_c2_23, x1=np.random.binomial(1,0.5), x2=np.random.binomial(1,0.5))
    gen[:,c3_23[0]], gen[:,c3_23[1]], x3_23_1,x3_23_2, rec3_23_1, rec3_23_2= make_baby(gen=gen, pos=pos, parent1=p2, parent2=p3,recomb_no1=no1_c3_23, recomb_interval1=interv1_c3_23, recomb_no2=no2_c3_23, recomb_interval2=interv2_c3_23, x1=np.random.binomial(1,0.5), x2=np.random.binomial(1,0.5))
    gen[:,c4_45[0]], gen[:,c4_45[1]], x4_45_1,x4_45_2,rec4_45_1, rec4_45_2= make_baby(gen=gen, pos=pos, parent1=p4, parent2=p5,recomb_no1=no1_c4_45, recomb_interval1=interv1_c4_45, recomb_no2=no2_c4_45, recomb_interval2=interv2_c4_45, x1=np.random.binomial(1,0.5), x2=np.random.binomial(1,0.5))
    gen[:,gc1_1011[0]], gen[:,gc1_1011[1]],x1c_1011_1,x1c_1011_2, regc1_1011_1, regc1_1011_2= make_baby(gen=gen, pos=pos, parent1=c3_23, parent2=c4_45,recomb_no1=no1_gc1_1011, recomb_interval1=interv1_gc1_1011, recomb_no2=no2_gc1_1011, recomb_interval2=interv2_gc1_1011, x1=np.random.binomial(1,0.5), x2=np.random.binomial(1,0.5))
    gen[:,ggc1_512[0]], gen[:,ggc1_512[1]],x1c_512_1,x1c_512_2, regc1_512_1, regc1_512_2= make_baby(gen=gen, pos=pos, parent1=gc1_1011, parent2=p6,recomb_no1=no1_ggc1_512, recomb_interval1=interv1_ggc1_512, recomb_no2=no2_ggc1_512, recomb_interval2=interv2_ggc1_512, x1=np.random.binomial(1,0.5), x2=np.random.binomial(1,0.5))
    gen[:,gggc1_613[0]], gen[:,gggc1_613[1]],x1c_613_1,x1c_613_2, regc1_613_1, regc1_613_2= make_baby(gen=gen, pos=pos, parent1=ggc1_512, parent2=p7,recomb_no1=no1_gggc1_613, recomb_interval1=interv1_gggc1_613, recomb_no2=no2_gggc1_613, recomb_interval2=interv2_gggc1_613, x1=np.random.binomial(1,0.5), x2=np.random.binomial(1,0.5))




    par=[p1,p2,p2,p3,p2,p3,p4,p5,c3_23,c4_45,gc1_1011,p6,ggc1_512,p7]
    chrm=list(range(unrelated,unrelated+related))
    rec=[rec1_12_1, rec1_12_2, rec2_23_1, rec2_23_2, rec3_23_1, rec3_23_2, rec4_45_1, rec4_45_2, regc1_1011_1, regc1_1011_2, regc1_512_1, regc1_512_2,regc1_613_1, regc1_613_2]
    x=[x1_12_1, x1_12_2, x2_23_1, x2_23_2, x3_23_1, x3_23_2, x4_45_1, x4_45_2, x1c_1011_1, x1c_1011_2, x1c_512_1, x1c_512_2, x1c_613_1, x1c_613_2]

    recomb_mat=pd.DataFrame(
            {'chromosome': chrm,
             'parents': par,
             'recomb_points': rec,
             'random_chr': x
             })

    path= 'sim_init' + '/'+ 'run' + str(RID) + '/' + chrm_no + 'sim_recomb'

    if os.path.exists(path):
        shutil.rmtree(path)
    os.mkdir(path)

    for c in range(len(chrm)):
        rfile_name='sim_init' + '/'+'run'+ str(RID)+'/'+ chrm_no +'sim_recomb/' + str(recomb_mat['chromosome'][c])+'par'+str(recomb_mat['parents'][c][0])+','+str(recomb_mat['parents'][c][1])+ '_' + str(recomb_mat['random_chr'][c])+'.csv'
        np.savetxt(fname=rfile_name, X=rec[c], delimiter=",")

    gen[:,(unrelated+related):(unrelated+related+identical)]=gen[:,0:identical]             #identical individuals

    #print(np.sum(gen,0))

    Dip=pd.DataFrame(
            data=gen,
            index=np.array(range(0,np.shape(gen)[0])),
            columns=np.array(range(0,np.shape(gen)[1]))
            )

    with pd.option_context('display.max_rows', len(Dip.index), 'display.max_columns', len(Dip.columns)):
                Dip.to_csv(diploid_sim, sep=',')


    np.savetxt(fname=out_hbdall, X=hbdall, delimiter=",")










def getAscertainment(ascertainment_scheme,dipfile,posfile,pos_out,chrm):
    dip=pd.read_csv(dipfile, sep=',',index_col=0,header=0)
    pos=np.loadtxt(posfile,dtype='float', delimiter = ",")

    dip.insert(loc=0,column='pos',value=pos)
    dip.insert(loc=0,column='chrm',value=chrm)

    x=int(ascertainment_scheme)*2

    toremove=['pos','chrm']
    toremove.extend(list(map(str, list(range(34,50)))))

    if x!=0:

        col=['chrm', 'pos']
        col.extend(list(range(0,x)))
        poly=pd.DataFrame(columns=col)

        poly['chrm']=dip['chrm']
        poly['pos']=dip['pos']

        poly.iloc[:,2:2+x]=dip.iloc[:,2+34:2+34+x].values       #ind 17,18,19 are high coverage, keep in mind that the first two columns are chrm and pos

        keepsites = pos[np.any(poly.iloc[:,2:] != poly.iloc[:,2][:, np.newaxis], 1)]



    elif x==0:
        df1=dip
        df=df1.drop(labels=toremove,axis=1)
        keepsites=pos

    np.savetxt(fname=pos_out, X=keepsites, delimiter=",")



def getReads(genf, coverage, iscontam, contaml, totalind, c_no, outf):

    gen=pd.read_csv(genf, sep=",", header=0,index_col=0)

    df= pd.DataFrame(columns=range(totalind*4))
    c_freq=np.sum(gen.loc[:,gen.columns[-int(c_no):]],axis=1)/int(c_no)
    coveragef=float(coverage)

    if int(iscontam)==0:
        contaml=np.zeros(len(contaml))


    for i in range(totalind):
        #ind1 = int(2*i)
        #ind2 = int(2*i+1)
        r1ea=int(4*i)
        r2ed=int(4*i + 1)
        r3ca=int(4*i + 2)
        r4cd=int(4*i + 3)
        #p= (gen.loc[:,str(ind1)] + gen.loc[:,str(ind2)])/2
        p=np.ones(len(gen)) * 0.5
        df[r1ea]=np.random.poisson(coveragef * (p) * (1-(contaml[i]/100)))
        df[r2ed]=np.random.poisson(coveragef * (p) * (1-(contaml[i]/100)))
        df[r3ca]=np.random.poisson(coveragef * (1-c_freq) * (contaml[i]/100))
        df[r4cd]=np.random.poisson(coveragef * (c_freq) * (contaml[i]/100))


    with pd.option_context('display.max_rows', len(df.index), 'display.max_columns', len(df.columns)):
             df.to_csv(outf, sep=',')


def finalReads(genf,read01,totalind, outf):

    gen=pd.read_csv(genf, sep=",", header=0,index_col=0)
    pre=pd.read_csv(read01, sep=",", header=0,index_col=0)

    for i in range(totalind):
            ind1 = int(2*i)
            ind2 = int(2*i+1)
            r1ea=int(4*i)
            r2ed=int(4*i + 1)

            true_ed=gen.loc[:,str(ind1)] * pre.loc[:,str(r1ea)] + gen.loc[:,str(ind2)] * pre.loc[:,str(r2ed)]
            true_ea=pre.loc[:,str(r1ea)] + pre.loc[:,str(r2ed)] - true_ed

            pre[str(r1ea)] = true_ea
            pre[str(r2ed)] = true_ed

    with pd.option_context('display.max_rows', len(pre.index), 'display.max_columns', len(pre.columns)):
                pre.to_csv(outf, sep=',')




def probDiffs(readsf, pos_ascf, pos_allf, chrm, fout_prob, fout_idDiff, fout_pshap, fout_idhap):

    reads=pd.read_csv(readsf, sep=',',index_col=0,header=0)
    pos_asc=np.loadtxt(pos_ascf,dtype='float', delimiter = ",")
    pos_all=np.loadtxt(pos_allf,dtype='float', delimiter = ",")

    l=int(np.shape(reads)[1])

    reads['chrm']=chrm
    reads['pos']=pos_all

    xx=reads['pos'].isin(pos_asc)

    reads=reads.loc[xx,:]
    reads.reset_index(drop=True, inplace=True)


    cols=list(map(str, list(range(l))))
    r1=reads[cols].to_numpy()
    r3=np.array(np.hsplit(r1,int(np.shape(r1)[1]/4)))


    p_ae=r3[:,:,0]
    p_de=r3[:,:,1]
    p_ac=r3[:,:,2]
    p_dc=r3[:,:,3]


    der= p_de + p_dc
    anc= p_ae + p_ac

    total=der+anc

    prob_ar= der / total


    idDiff=scp.comb(der,1) * scp.comb(anc,1)/ scp.comb(total, 2)



    pshap=np.random.binomial(1,prob_ar)
    pshap[pshap<0]=-9

    idhap=np.random.binomial(1, idDiff)
    idhap[total<2]=-9

    prob_ar[total==0]=-9
    idDiff[total<2]=-9

    np.savetxt(fname=fout_prob, X=prob_ar.T, delimiter=",")
    np.savetxt(fname=fout_idDiff, X=idDiff.T, delimiter=",")
    np.savetxt(fname=fout_pshap, X=pshap.T, delimiter=",")
    np.savetxt(fname=fout_idhap, X=idhap.T, delimiter=",")



def getPseudohap(genotype, genotype_i, outfile, outfile_i, ind, chrm, position, position_i):
    gen= pd.read_csv(genotype, sep=",", index_col=0, header=0)
    #print(gen)
    pos= np.loadtxt(position,dtype='float', delimiter = ",")

    diploid= np.zeros((int(np.shape(gen)[1]/2),2))      #making array of chromosomes that belong to same individual
    for i in range(0,np.shape(gen)[1],2):
        diploid[int(i/2),:]=(i,i+1)

    i=int(ind)

    ind1=str((diploid[i]).astype(int)[0])
    ind2=str((diploid[i]).astype(int)[1])

    ps12=gen[[str(ind1),str(ind2)]].to_numpy()

    a=np.zeros(np.shape(gen)[0])

    ran=np.random.choice([0, 1], size=(np.shape(gen)[0]), p=[0.5, 0.5])
    for po in range(np.shape(gen)[0]):
        a[po]= ps12[po,ran[po]]


    df = pd.DataFrame(
        {'allele': a,
         'chrm': chrm,
         'pos': pos
         })

    df['allele']=df['allele'].replace(1,'G')
    df['allele']=df['allele'].replace(0,'A')

    with pd.option_context('display.max_rows', len(df.index), 'display.max_columns', len(df.columns)):
        df.to_csv(outfile, sep=',')


    gen_i= pd.read_csv(genotype_i, sep=",", index_col=0, header=0)
    #print(gen_i)
    pos_i= np.loadtxt(position_i,dtype='float', delimiter = ",")


    ps12_i=gen_i[[str(ind1),str(ind2)]].to_numpy()

    a_i=np.zeros(np.shape(gen_i)[0])

    for po in range(np.shape(gen_i)[0]):
        a_i[po]= ps12_i[po,ran[po]]

    df_i = pd.DataFrame(
        {'allele': a_i,
         'chrm': chrm,
         'pos': pos_i
         })

    df_i['allele']=df_i['allele'].replace(1,'G')
    df_i['allele']=df_i['allele'].replace(0,'A')

    with pd.option_context('display.max_rows', len(df_i.index), 'display.max_columns', len(df_i.columns)):
        df_i.to_csv(outfile_i, sep=',')




def get_transition(gstate, sibstate, sibAfile, grandAfile, e, no_st):


    states=np.loadtxt(fname=sibstate, dtype='float', delimiter=',')
    stateg=np.loadtxt(fname=gstate, dtype='float', delimiter=',')
    state_array=[states,stateg]
    out=[sibAfile,grandAfile]

    for rel in range(len(state_array)):
        state=state_array[rel]
        sibA=np.zeros((no_st,no_st))

        for i in range(len(state)-1):
            if state[i]==0.75:
                if state[i+1]==0.75:
                    sibA[0,0]=sibA[0,0]+1
                elif state[i+1]==1:
                    sibA[0,1]=sibA[0,1]+1
                elif state[i+1]==0.5:
                    sibA[0,2]=sibA[0,2]+1
            elif state[i]==1:
                if state[i+1]==0.75:
                    sibA[1,0]=sibA[1,0]+1
                elif state[i+1]==1:
                    sibA[1,1]=sibA[1,1]+1
                elif state[i+1]==0.5:
                    sibA[1,2]=sibA[1,2]+1
            if state[i]==0.5:
                if state[i+1]==0.75:
                    sibA[2,0]=sibA[2,0]+1
                elif state[i+1]==1:
                    sibA[2,1]=sibA[2,1]+1
                elif state[i+1]==0.5:
                    sibA[2,2]=sibA[2,2]+1

        sibA=sibA/sibA.sum(axis=1)[:,None]

        if rel==1:              #for grand parent child pair
            sibA[2,:]=[0.5-e,0.5-e,2*e]
            sibA[:2,2]=sibA[:2,2]+(2*e)
            sibA[0,:2]=sibA[0,:2]-e
            sibA[1,:2]=sibA[1,:2]-e
        #print(state_array[rel])
        df_A=pd.DataFrame(sibA)

        #print(df_A)

        with pd.option_context('display.max_rows', len(df_A.index), 'display.max_columns', len(df_A.columns)):
                    df_A.to_csv(out[rel], sep=',', header=False, index=False)



def avgA(siblistA, grandlistA, e, l, sibAfinal,grandAfinal,pcAfinal,idAfinal,unAfinal):



    Ag=np.zeros((l,l))
    Asib=np.zeros((l,l))

    for i in range(len(siblistA)):
        gmat=pd.read_csv(grandlistA[i], sep=",", header=None,index_col=False)
        sibmat=pd.read_csv(siblistA[i], sep=",", header=None,index_col=False)
        Ag=Ag+gmat
        Asib=Asib+sibmat

    Ag=pd.DataFrame(Ag/(i+1))
    Asib=pd.DataFrame(Asib/(i+1))
    Apc=pd.DataFrame(np.array([[1-(2*e),e,e],[1-(2*e),e,e],[1-(2*e),e,e]]))
    Aun=pd.DataFrame(np.array([[e,1-(2*e),e],[e,1-(2*e),e],[e,1-(2*e),e]]))
    Aid=pd.DataFrame(np.array([[e,e,1-(2*e)],[e,e,1-(2*e)],[e,e,1-(2*e)]]))

    with pd.option_context('display.max_rows', len(Ag.index), 'display.max_columns', len(Ag.columns)):
                Ag.to_csv(grandAfinal, sep=',', header=False, index=False)

    with pd.option_context('display.max_rows', len(Asib.index), 'display.max_columns', len(Asib.columns)):
                Asib.to_csv(sibAfinal, sep=',', header=False, index=False)

    with pd.option_context('display.max_rows', len(Apc.index), 'display.max_columns', len(Apc.columns)):
                Apc.to_csv(pcAfinal, sep=',', header=False, index=False)

    with pd.option_context('display.max_rows', len(Aun.index), 'display.max_columns', len(Aun.columns)):
                        Aun.to_csv(unAfinal, sep=',', header=False, index=False)

    with pd.option_context('display.max_rows', len(Aid.index), 'display.max_columns', len(Aid.columns)):
                        Aid.to_csv(idAfinal, sep=',', header=False, index=False)
