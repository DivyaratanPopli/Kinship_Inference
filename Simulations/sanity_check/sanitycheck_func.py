import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import figure
import os
import math

#preparing the recombination points file
def prep_rec(rdir):

    files=np.sort(os.listdir(rdir))
    l=len(files)
    chrm=np.zeros(l)
    x=np.zeros(l)
    par=[]
    rec=[]
    for f in range(l):
        rname=rdir+str(files[f])
        rfile=np.loadtxt(rname,dtype='float', delimiter = ",")
        chrm[f]=rname.split('par')[0].split('recomb/')[1]

        x[f]=rname.split('.csv')[0].split('_')[3]
        par.append(tuple(rname.split('par')[1].split('.csv')[0].split('_')[0].split(',')))
        rec.append(rfile)
    chrm=chrm.astype(int)
    x=x.astype(int)
    recomb_mat=pd.DataFrame(
                {'chromosome': chrm,
                 'parents': par,
                 'recomb_points': rec,
                 'random_chrm': x
                 })
    return recomb_mat


#################################################### Finding recombination points

################################ grandparent-grandchild

def plot_gp(length, interval,filenamed,filenamet,rdir,g1,g2,p1,p2,outimg,outstate,rel23, ninds):

    g1=str(g1)
    g2=str(g2)
    p1=str(p1)
    p2=str(p2)

    win=math.ceil(length/interval)

    rec=prep_rec(rdir)

    j=4
    i=0
    state_32=np.zeros(int(rec.loc[j,"recomb_points"][-1]))

    while i<np.shape(rec.loc[j,"recomb_points"])[0]-1:

        if i%2==0:
            state_32[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=rec.loc[j,"parents"][rec['random_chrm'][j]]
        elif i%2==1:
            state_32[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=rec.loc[j,"parents"][abs(rec['random_chrm'][j]-1)]
        i=i+1


    #plt.plot(state_32)

    j=5
    i=0
    state_33=np.zeros(int(rec.loc[j,"recomb_points"][-1]))

    while i<np.shape(rec.loc[j,"recomb_points"])[0]-1:

        if i%2==0:
            state_33[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=rec.loc[j,"parents"][rec['random_chrm'][j]]
        elif i%2==1:
            state_33[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=rec.loc[j,"parents"][abs(rec['random_chrm'][j]-1)]
        i=i+1


    #plt.plot(state_33)

    j=8
    i=0
    state_36=np.zeros(int(rec.loc[j,"recomb_points"][-1]))
    state_array=np.vstack((state_32,state_33))

    while i<np.shape(rec.loc[j,"recomb_points"])[0]-1:

        if i%2==0:
            state_36[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=state_array[rec['random_chrm'][j],int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]    # for the old method without randomisation, see earlier version
        elif i%2==1:
            state_36[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=state_array[abs(rec['random_chrm'][j]-1), int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]
        i=i+1


    if rel23==3 or rel23==4 or rel23==5:
        j=6
        i=0
        state_34=np.zeros(int(rec.loc[j,"recomb_points"][-1]))

        while i<np.shape(rec.loc[j,"recomb_points"])[0]-1:

            if i%2==0:
                state_34[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=rec.loc[j,"parents"][rec['random_chrm'][j]]
            elif i%2==1:
                state_34[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=rec.loc[j,"parents"][abs(rec['random_chrm'][j]-1)]
            i=i+1


        #plt.plot(state_32)

        j=7
        i=0
        state_35=np.zeros(int(rec.loc[j,"recomb_points"][-1]))

        while i<np.shape(rec.loc[j,"recomb_points"])[0]-1:

            if i%2==0:
                state_35[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=rec.loc[j,"parents"][rec['random_chrm'][j]]
            elif i%2==1:
                state_35[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=rec.loc[j,"parents"][abs(rec['random_chrm'][j]-1)]
            i=i+1


        j=9
        i=0
        state_37=np.zeros(int(rec.loc[j,"recomb_points"][-1]))
        state_array=np.vstack((state_34,state_35))

        while i<np.shape(rec.loc[j,"recomb_points"])[0]-1:

            if i%2==0:
                state_37[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=state_array[rec['random_chrm'][j],int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]    # for the old method without randomisation, see earlier version
            elif i%2==1:
                state_37[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=state_array[abs(rec['random_chrm'][j]-1), int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]
            i=i+1


        j=10
        i=0
        state_38=np.zeros(int(rec.loc[j,"recomb_points"][-1]))
        state_array=np.vstack((state_36,state_37))

        while i<np.shape(rec.loc[j,"recomb_points"])[0]-1:

            if i%2==0:
                state_38[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=state_array[rec['random_chrm'][j],int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]    # for the old method without randomisation, see earlier version
            elif i%2==1:
                state_38[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=state_array[abs(rec['random_chrm'][j]-1), int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]
            i=i+1

        if rel23==4 or rel23==5:

            state_39=np.ones(int(rec.loc[j,"recomb_points"][-1]))*10

            j=12
            i=0
            state_40=np.zeros(int(rec.loc[j,"recomb_points"][-1]))
            state_array=np.vstack((state_38,state_39))

            while i<np.shape(rec.loc[j,"recomb_points"])[0]-1:

                if i%2==0:
                    state_40[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=state_array[rec['random_chrm'][j],int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]    # for the old method without randomisation, see earlier version
                elif i%2==1:
                    state_40[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=state_array[abs(rec['random_chrm'][j]-1), int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]
                i=i+1

        if rel23==5:

            j=1
            i=0
            state_41=np.zeros(int(rec.loc[j,"recomb_points"][-1]))

            while i<np.shape(rec.loc[j,"recomb_points"])[0]-1:

                if i%2==0:
                    state_41[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=rec.loc[j,"parents"][rec['random_chrm'][j]]
                elif i%2==1:
                    state_41[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=rec.loc[j,"parents"][abs(rec['random_chrm'][j]-1)]
                i=i+1




    if g1=='1':
        state_par=np.copy(state_32)
    elif g1=='2':
        state_par=np.copy(state_33)
    elif g1=='8':
        state_par=np.copy(state_41)

    if rel23==2:
        wh_state=np.copy(state_36)
    elif rel23==3:
        wh_state=np.copy(state_38)
    elif rel23==4 or rel23==5:
        wh_state=np.copy(state_40)


    state2or3=np.ones(len(wh_state))*0.1
    state2or3[wh_state==state_par]=0
    state2or3[wh_state!=state_par]=1



    ind=(np.where(state2or3[:-1] != state2or3[1:])[0]*win/length).astype(int)
    index=np.zeros(np.shape(ind)[0]+2)
    index[0]=0
    index[-1]=win
    index[1:-1]=ind
    index=index.astype(int)

    diff= np.loadtxt(filenamed,dtype='float', delimiter = ",")
    total=np.loadtxt(filenamet,dtype='float', delimiter = ",")

    xxx=[]
    for x in range(ninds-1):
        for y in range(x+1,ninds):
            xxx.append(str(x)+'_'+str(y))
    xx=np.array(xxx)
    ig=np.where(xx==str(g1)+'_'+str(g2))[0][0]
    ip=np.where(xx==str(p1)+'_'+str(p2))[0][0]

    print("ig", ig, "ip", ip, "gpair=",str(g1)+'_'+str(g2),"ppair=",str(p1)+'_'+str(p2))
    data=diff[:,ig]/total[:,ig]
    data1=diff[:,ip]/total[:,ip]


    #plt.figure(figsize=(10,10))
    xax=np.array(range(win))

    state_final= state2or3[(xax/win*length).astype(int)+int(interval/2)]
    state_final[state_final==0]=3/4

    #np.savetxt(fname=outstate, X=state_final, delimiter=',')
    np.savetxt(fname=outstate, X=state_final, delimiter=',')

    fig, ax = plt.subplots(nrows=2, ncols=1)
    ax[0].plot(data,'blue')
    ax[0].plot(data1,'purple')

    ax[1].plot(xax,state_final, color='black', linestyle='dashed',linewidth=1.5, markersize=7)


    plt.savefig(outimg)

    plt.close()



    return index, state_final


#######################################################################siblings

#making sibling 1 as fragments of the parents
#length=5e7


def plot_sib(length, interval,filenamed,filenamet,rdir,p11,p12,p21,p22, outimg,outstate, ninds, rel):


    win=math.ceil(length/interval)

    rec=prep_rec(rdir)

    j=p11             #row number for rec (index for the chromosome)
    i=0             #index for the shared fragment
    state_30=np.zeros(int(rec.loc[j,"recomb_points"][-1]))

    while i<np.shape(rec.loc[j,"recomb_points"])[0]-1:

        if i%2==0:
            state_30[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=rec.loc[j,"parents"][rec['random_chrm'][j]]
        elif i%2==1:
            state_30[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=rec.loc[j,"parents"][abs(rec['random_chrm'][j]-1)]
        i=i+1


    #plt.plot(state_32)

    j=p12
    i=0
    state_31=np.zeros(int(rec.loc[j,"recomb_points"][-1]))

    while i<np.shape(rec.loc[j,"recomb_points"])[0]-1:

        if i%2==0:
            state_31[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=rec.loc[j,"parents"][rec['random_chrm'][j]]
        elif i%2==1:
            state_31[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=rec.loc[j,"parents"][abs(rec['random_chrm'][j]-1)]
        i=i+1


    #making sibling 2 as fragments of parents



    j=p21             #row number for rec (index for the chromosome)
    i=0             #index for the shared fragment
    state_32=np.zeros(int(rec.loc[j,"recomb_points"][-1]))

    while i<np.shape(rec.loc[j,"recomb_points"])[0]-1:

        if i%2==0:
            state_32[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=rec.loc[j,"parents"][rec['random_chrm'][j]]
        elif i%2==1:
            state_32[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=rec.loc[j,"parents"][abs(rec['random_chrm'][j]-1)]
        i=i+1


    #plt.plot(state_32)

    j=p22
    i=0
    state_33=np.zeros(int(rec.loc[j,"recomb_points"][-1]))

    while i<np.shape(rec.loc[j,"recomb_points"])[0]-1:

        if i%2==0:
            state_33[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=rec.loc[j,"parents"][rec['random_chrm'][j]]
        elif i%2==1:
            state_33[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=rec.loc[j,"parents"][abs(rec['random_chrm'][j]-1)]
        i=i+1


    if rel=='avu':

        j=8
        i=0
        state_40=np.zeros(int(rec.loc[j,"recomb_points"][-1]))
        state_array=np.vstack((state_32,state_33))

        while i<np.shape(rec.loc[j,"recomb_points"])[0]-1:

            if i%2==0:
                state_40[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=state_array[rec['random_chrm'][j],int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]    # for the old method without randomisation, see earlier version
            elif i%2==1:
                state_40[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=state_array[abs(rec['random_chrm'][j]-1), int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]
            i=i+1



    if rel=='sib':
        sib_state=np.zeros(np.shape(state_32)[0])

        sib_state[state_32!=state_30]=sib_state[state_32!=state_30]+1
        sib_state[state_32!=state_31]=sib_state[state_32!=state_31]+1
        sib_state[state_33!=state_30]=sib_state[state_33!=state_30]+1
        sib_state[state_33!=state_31]=sib_state[state_33!=state_31]+1

        sib_state=sib_state/4
        pair='9_10'

    elif rel=='avu':
        sib_state=np.zeros(np.shape(state_40)[0])
        sib_state[(state_40==state_30) + (state_40==state_31)]=0.75
        sib_state[~((state_40==state_30) + (state_40==state_31))]=1
        pair='9_12'

    pairp='1_9'

    ind=(np.where(sib_state[:-1] != sib_state[1:])[0]*win/length).astype(int)
    index=np.zeros(np.shape(ind)[0]+2)
    index[0]=0
    index[-1]=win
    index[1:-1]=ind
    index=index.astype(int)


    diff= np.loadtxt(filenamed,dtype='float', delimiter = ",")
    total=np.loadtxt(filenamet,dtype='float', delimiter = ",")

    xxx=[]
    for x in range(ninds-1):
        for y in range(x+1,ninds):
            xxx.append(str(x)+'_'+str(y))


    xx=np.array(xxx)
    isa=np.where(xx==pair)[0][0]
    ip=np.where(xx==pairp)[0][0]

    print("isa", isa, "ip", ip, "spair=",pair,"ppair=",pairp)
    data=diff[:,isa]/total[:,isa]
    data1=diff[:,ip]/total[:,ip]

    xax=np.array(range(win))

    final_state=sib_state[(xax/win*length).astype(int)+int(interval/2)]
    np.savetxt(fname=outstate, X=final_state, delimiter=',')

    fig, ax = plt.subplots(nrows=2, ncols=1)
    ax[0].plot(data, 'blue')
    ax[0].plot(data1,'purple')

    ax[1].plot(xax,final_state, color='black', linestyle='dashed',linewidth=1.5, markersize=7)

    plt.savefig(outimg)
    plt.close()



    return index, final_state



def plot_halfsib(length, interval,filenamed, filenamet, rdir, outimg,outstate, ninds):

    win=math.ceil(length/interval)

    rec=prep_rec(rdir)

    j=1             #row number for rec (index for the chromosome)
    i=0             #index for the shared fragment
    state_31=np.zeros(int(rec.loc[j,"recomb_points"][-1]))

    while i<np.shape(rec.loc[j,"recomb_points"])[0]-1:

        if i%2==0:
            state_31[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=rec.loc[j,"parents"][rec['random_chrm'][j]]
        elif i%2==1:
            state_31[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=rec.loc[j,"parents"][abs(rec['random_chrm'][j]-1)]
        i=i+1


    j=2             #row number for rec (index for the chromosome)
    i=0             #index for the shared fragment
    state_32=np.zeros(int(rec.loc[j,"recomb_points"][-1]))

    while i<np.shape(rec.loc[j,"recomb_points"])[0]-1:

        if i%2==0:
            state_32[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=rec.loc[j,"parents"][rec['random_chrm'][j]]
        elif i%2==1:
            state_32[int(rec.loc[j,"recomb_points"][i]):int(rec.loc[j,"recomb_points"][i+1])]=rec.loc[j,"parents"][abs(rec['random_chrm'][j]-1)]
        i=i+1

    state_hs=np.zeros(int(rec.loc[j,"recomb_points"][-1]))

    for fi in range(len(state_hs)):
        if state_31[fi]!=state_32[fi]:
            state_hs[fi]=1


    ind=(np.where(state_hs[:-1] != state_hs[1:])[0]*win/length).astype(int)
    index=np.zeros(np.shape(ind)[0]+2)
    index[0]=0
    index[-1]=win
    index[1:-1]=ind
    index=index.astype(int)

    diff= np.loadtxt(filenamed,dtype='float', delimiter = ",")
    total=np.loadtxt(filenamet,dtype='float', delimiter = ",")

    pair='8_9'
    pairp='1_9'
    xxx=[]
    for x in range(ninds-1):
        for y in range(x+1,ninds):
            xxx.append(str(x)+'_'+str(y))


    xx=np.array(xxx)
    ihs=np.where(xx==pair)[0][0]
    ip=np.where(xx==pairp)[0][0]

    print("ihs", ihs, "ip", ip, "hspair=",pair,"ppair=",pairp)
    data=diff[:,ihs]/total[:,ihs]
    data1=diff[:,ip]/total[:,ip]


    xax=np.array(range(win))

    state_final= state_hs[(xax/win*length).astype(int)+int(interval/2)]
    state_final[state_final==0]=3/4

    np.savetxt(fname=outstate, X=state_final, delimiter=',')

    fig, ax = plt.subplots(nrows=2, ncols=1)
    ax[0].plot(data, 'blue')
    ax[0].plot(data1, 'purple')

    ax[1].plot(xax,state_final, color='black', linestyle='dashed',linewidth=1.5, markersize=7)


    plt.savefig(outimg)

    plt.close()


    return index, state_final


def merge_all(infiles,outfile):

    #merging the data files
    columns=['chrom', 'dis', 'count']
    mer=pd.DataFrame(columns=columns)



    for fi in range(len(infiles)):

        temp=pd.read_csv(infiles[fi], sep=",", index_col=0,header=0)
        if infiles[fi].split('chrm')[1].split('/')[0] != str(fi+1):
            print("files are not sorted for merge_all fn")
        ch_array=(np.ones(len(temp))*(fi+1)).astype(int)
        temp.insert(loc=0,column='chrom', value=ch_array)
        mer=mer.append(temp)



    index=np.array(range(len(mer)))
    mer['index']=index
    mer=mer.set_index(keys=index)
    mer=mer.drop('index', axis=1)



    with pd.option_context('display.max_rows', len(mer.index), 'display.max_columns', len(mer.columns)):
                mer.to_csv(outfile, sep=',')




def make_plots(grandstate,sibstate,hsibstate,deg3state,avustate,deg4state,deg5state,diff,total,fig_g,fig_s,fig_h,fig_d3,fig_a,fig_d4,fig_d5, ninds):

    obsd= np.loadtxt(diff,dtype='float', delimiter = ",")
    obst= np.loadtxt(total,dtype='float', delimiter = ",")

    xxx=[]
    for x in range(ninds-1):
        for y in range(x+1,ninds):
            xxx.append(str(x)+'_'+str(y))

    xx=np.array(xxx)

    inds,indg,indp,indh,indd3,inda,indd4,indd5=np.where(xx=='9_10')[0][0],np.where(xx=='1_12')[0][0],np.where(xx=='1_10')[0][0],np.where(xx=='8_9')[0][0],np.where(xx=='1_13')[0][0],np.where(xx=='9_12')[0][0],np.where(xx=='1_14')[0][0],np.where(xx=='8_14')[0][0]
    datas=obsd[:,inds]/obst[:,inds]
    datag=obsd[:,indg]/obst[:,indg]
    datap=obsd[:,indp]/obst[:,indp]
    datah=obsd[:,indh]/obst[:,indh]
    datad3=obsd[:,indd3]/obst[:,indd3]
    dataa=obsd[:,inda]/obst[:,inda]
    datad4=obsd[:,indd4]/obst[:,indd4]
    datad5=obsd[:,indd5]/obst[:,indd5]

    if len(datas)==len(datap) and len(datag)==len(datap) and len(datag)==len(datad3) and len(datag)==len(datah) and len(datag)==len(dataa) and len(datag)==len(datad4) and len(datag)==len(datad5):
        win=len(datas)
    else:
        print("Warning:there is some problem with missing windows")
    xax=np.array(range(win))

    g_state=np.loadtxt(fname=grandstate, dtype='float', delimiter=',')
    s_state=np.loadtxt(fname=sibstate, dtype='float', delimiter=',')
    h_state=np.loadtxt(fname=hsibstate, dtype='float', delimiter=',')
    d3_state=np.loadtxt(fname=deg3state, dtype='float', delimiter=',')
    a_state=np.loadtxt(fname=avustate, dtype='float', delimiter=',')
    d4_state=np.loadtxt(fname=deg4state, dtype='float', delimiter=',')
    d5_state=np.loadtxt(fname=deg5state, dtype='float', delimiter=',')

    print('xax',xax)
    print('gstate',g_state)
    #plotting grand state
    alldata=[datag,datas,datah,datad3,dataa,datad4,datad5]
    allst=[g_state,s_state,h_state,d3_state,a_state,d4_state,d5_state]
    fig_name=[fig_g,fig_s,fig_h,fig_d3,fig_a,fig_d4,fig_d5]

    for st in range(len(alldata)):
        fig, ax = plt.subplots(nrows=2, ncols=1)
        ax[0].plot(alldata[st], 'blue')
        if st==5:
            datap=datad3.copy()
        ax[0].plot(datap,'purple')

        ax[1].plot(xax,allst[st], color='black', linestyle='dashed',linewidth=1.5, markersize=7)

        plt.savefig(fig_name[st])
        plt.close()
