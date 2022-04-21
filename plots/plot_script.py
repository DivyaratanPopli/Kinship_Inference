import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def ROH_plotf(fold,outfold):

    a=0
    c=0
    r=1
    i=1
    ind=0
    cov=4
    f=0
    rohf=fold+"hbd_hmm/contam%s/inbred%s/run%s/coverage%s/filtered%s/asc%s_res/gamma_hapProbs/pw_%s.csv.gz" %(c,i,r,cov,f,a,ind)
    roh=pd.read_csv(rohf, sep=",", header=0,index_col=0)
    roh.loc[roh['g_noin']>1,'g_noin']=np.nan

    difff=fold+"contam%s/inbred%s/run%s/coverage%s/asc%s/inputMode_idhapProbs_fil%s/id_wind.csv.gz" %(c,i,r,cov,a,f)
    totalf=fold+"contam%s/inbred%s/run%s/coverage%s/asc%s/inputMode_idhapProbs_fil%s/id_wint.csv.gz" %(c,i,r,cov,a,f)

    diff=np.loadtxt(difff, dtype='float', delimiter = ",")[:,ind]
    total=np.loadtxt(totalf, dtype='float', delimiter = ",")[:,ind]


    truedf=fold+"post_sim/inbred%s/run%s/hbd_merged/diff_all.csv.gz" %(i,r)
    truetf=fold+"post_sim/inbred%s/run%s/hbd_merged/total_all.csv.gz" %(i,r)

    hbdd=np.loadtxt(truedf, dtype='float', delimiter = ",")[:,ind]
    hbdt=np.loadtxt(truetf, dtype='float', delimiter = ",")[:,ind]

    #plt.plot(diff/total)
    #plt.plot(roh['g_noin']/20)
    #plt.plot(hbdd/hbdt/20)

    df1=pd.DataFrame(
        columns=['Chrm', 'Prop', 'Coverage',
           'True_ROH', 'Model_pred', 'acr','Wins']
        )


    for pl in range(4):
        if pl==0:
            a=0
            cnt=0
            inb=0
        elif pl==1:
            a=0
            cnt=0
            inb=1
        elif pl==2:
            a=0
            cnt=1
            inb=1
        elif pl==3:
            a=2
            cnt=0
            inb=1

        acr=str(a)+ '_'+ str(cnt)+ '_'+ str(inb)
        for cov in [4,0.2,0.1]:

            rohf=fold+"hbd_hmm/contam%s/inbred%s/run%s/coverage%s/filtered%s/asc%s_res/gamma_hapProbs/pw_%s.csv.gz" %(cnt,inb,r,cov,f,a,ind)
            roh=pd.read_csv(rohf, sep=",", header=0,index_col=0)
            roh.loc[roh['g_noin']>1,'g_noin']=np.nan

            difff=fold+"contam%s/inbred%s/run%s/coverage%s/asc%s/inputMode_idhapProbs_fil%s/id_wind.csv.gz" %(cnt,inb,r,cov,a,f)
            totalf=fold+"contam%s/inbred%s/run%s/coverage%s/asc%s/inputMode_idhapProbs_fil%s/id_wint.csv.gz" %(cnt,inb,r,cov,a,f)

            diff=np.loadtxt(difff, dtype='float', delimiter = ",")[:,ind]
            total=np.loadtxt(totalf, dtype='float', delimiter = ",")[:,ind]


            truedf=fold+"post_sim/inbred%s/run%s/hbd_merged/diff_all.csv.gz" %(inb,r)
            truetf=fold+"post_sim/inbred%s/run%s/hbd_merged/total_all.csv.gz" %(inb,r)

            hbdd=np.loadtxt(truedf, dtype='float', delimiter = ",")[:,ind]
            hbdt=np.loadtxt(truetf, dtype='float', delimiter = ",")[:,ind]

            df=pd.DataFrame(
                    columns=['Chrm', 'Prop', 'Coverage',
                         'True_ROH', 'Model_pred', 'acr','Wins']
                            )
            df[['Chrm','Model_pred']]=roh[['chrom','g_noin']]
            df['Prop']=diff/total
            df['Coverage']=cov
            df['True_ROH']=hbdd/hbdt
            df['acr']=acr
            df['Wins']=np.array(list(range(len(roh))))+1


            df1=df1.append(df)

    with pd.option_context('display.max_rows', len(df1.index), 'display.max_columns', len(df1.columns)):
                df1.to_csv(outfold+'plot_ROH.csv', sep=',')


def IBD_plotf(fold,outfold):
    difff = fold+"contam0/inbred0/run1/coverage4/asc0/inputMode_hapProbs_fil0/remHighdiv_diff.csv.gz"
    totalf = fold+"contam0/inbred0/run1/coverage4/asc0/inputMode_hapProbs_fil0/remHighdiv_total.csv.gz"

    inds=list(range(17))
    listf=[]
    for i, l1 in enumerate(inds):
        for l2 in inds[(i+1):]:
            s = '%s_%s' %(str(l1).strip("['']"), str(l2).strip("['']"))
            listf.append(s)

    listf=np.array(listf)
    pair=np.where(listf=='9_10')[0][0]

    diff=np.loadtxt(difff, dtype='float', delimiter = ",")[:,pair]
    total=np.loadtxt(totalf, dtype='float', delimiter = ",")[:,pair]
    plt.plot(diff/total)

    res_unf=fold+"hmm_numba_corrA_relSim_fast_theoretical_sib_grA/contam0/inbred0/run1/coverage4/filtered0/asc0/relation_un_file/res_inphapProbs/pw_9_10.csv.gz"
    res_deg3f=fold+"hmm_numba_corrA_relSim_fast_theoretical_sib_grA/contam0/inbred0/run1/coverage4/filtered0/asc0/relation_deg3_file/res_inphapProbs/pw_9_10.csv.gz"
    res_secf=fold+"hmm_numba_corrA_relSim_fast_theoretical_sib_grA/contam0/inbred0/run1/coverage4/filtered0/asc0/relation_gr_file/res_inphapProbs/pw_9_10.csv.gz"
    res_pcf=fold+"hmm_numba_corrA_relSim_fast_theoretical_sib_grA/contam0/inbred0/run1/coverage4/filtered0/asc0/relation_pc_file/res_inphapProbs/pw_9_10.csv.gz"
    res_sibf=fold+"hmm_numba_corrA_relSim_fast_theoretical_sib_grA/contam0/inbred0/run1/coverage4/filtered0/asc0/relation_sib_file/res_inphapProbs/pw_9_10.csv.gz"
    res_idf=fold+"hmm_numba_corrA_relSim_fast_theoretical_sib_grA/contam0/inbred0/run1/coverage4/filtered0/asc0/relation_id_file/res_inphapProbs/pw_9_10.csv.gz"

    res_un=np.loadtxt(res_unf, dtype='float', delimiter = ",")
    res_deg3=np.loadtxt(res_deg3f, dtype='float', delimiter = ",")
    res_sec=np.loadtxt(res_secf, dtype='float', delimiter = ",")
    res_pc=np.loadtxt(res_pcf, dtype='float', delimiter = ",")
    res_sib=np.loadtxt(res_sibf, dtype='float', delimiter = ",")
    res_id=np.loadtxt(res_idf, dtype='float', delimiter = ",")


    trueibdf=fold+"sanity_check/inbred0/run1/coverage4/contam0/asc0/inputMode_hapProbs_fil0/sib_state_all.gz"
    trueibd=np.loadtxt(trueibdf, dtype='float', delimiter = ",")

    df_ibd=pd.DataFrame(
        {
         'Wins':np.array(list(range(220)))+1,
         'Prop':diff/total,
         'True':trueibd,
         'Unrelated':res_un,
         'degree3':res_deg3,
         'Degree2':res_sec,
         'Parent-Child':res_pc,
         'Siblings':res_sib,
         'Identical':res_id
         })

    with pd.option_context('display.max_rows', len(df_ibd.index), 'display.max_columns', len(df_ibd.columns)):
                df_ibd.to_csv(outfold+'plot_IBD.csv', sep=',')



def read_manipulation(t1):
    t=t1.copy()
    t=t.rename(columns = {'gr':'sec','un':'un_all'})
    t.loc['un_all',:]=(t.loc['un_all',:]* 94 + t.loc['deg3',:]* 9)/103
    t=t.drop('deg3',axis=0)
    t = t.reindex(index=['un_all','sec', 'fir', 'id'], columns=['un_all', 'sec', 'fir', 'id', 'NF'])


    t.loc['un_all',:]=t.loc['un_all',:] * 103
    t.loc['sec',:]=t.loc['sec',:] * 12
    t.loc['fir',:]=t.loc['fir',:] * 19
    t.loc['id',:]=t.loc['id',:] * 2
    t=t.round(decimals=3)

    return t



def kin_manipulation(t1):

    t=t1.copy()
    t['un_all']=t['un']+t['deg5']+t['deg4']+t['deg3']
    t=t.drop(['un','deg5','deg4','deg3'],axis=1)

    t.loc['un_all',:]=(t.loc['un_all',:]* 94 + t.loc['deg3',:]* 9)/103
    t=t.drop('deg3',axis=0)

    t['fir']=t['sib']+t['pc']
    t=t.drop(['pc','sib'],axis=1)

    t['sec']=t['gr']+t['avu']+t['hsib']
    t=t.drop(['gr','avu','hsib'],axis=1)


    t = t.reindex(index=['un_all', 'sec', 'fir', 'id'], columns=['un_all', 'sec', 'fir', 'id', 'NF'])

    t.loc['un_all',:]=t.loc['un_all',:] * 103
    t.loc['sec',:]=t.loc['sec',:] * 12
    t.loc['fir',:]=t.loc['fir',:] * 19
    t.loc['id',:]=t.loc['id',:] * 2

    return t

def kin_manipulation3(t1):

    t=t1.copy()
    t['un_all']=t['un']+t['deg5']+t['deg4']
    t=t.drop(['un','deg5','deg4'],axis=1)

    t['fir']=t['sib']+t['pc']
    t=t.drop(['pc','sib'],axis=1)

    t['sec']=t['gr']+t['avu']+t['hsib']
    t=t.drop(['gr','avu','hsib'],axis=1)


    t = t.reindex(index=['un_all', 'deg3', 'sec', 'fir', 'id'], columns=['un_all', 'deg3', 'sec', 'fir', 'id', 'NF'])

    t.loc['un_all',:]=t.loc['un_all',:] * 94
    t.loc['deg3',:]=t.loc['deg3',:] * 9
    t.loc['sec',:]=t.loc['sec',:] * 12
    t.loc['fir',:]=t.loc['fir',:] * 19
    t.loc['id',:]=t.loc['id',:] * 2

    return t

def all_manipulation(t1):
    t=t1.copy()
    t.loc['un_all',:]=(t.loc['un',:]* 87 + t.loc['deg5',:]* 1 + t.loc['deg4',:]* 6)/(94)
    t=t.drop(['un','deg4','deg5'],axis=0)
    t.loc['sec',:]=(t.loc['gr',:]* 9 + t.loc['hsib',:]* 2 + t.loc['avu',:])/12
    t=t.drop(['gr','hsib','avu'],axis=0)
    t.loc['fir',:]=(t.loc['pc',:]* 18 + t.loc['sib',:]* 1)/19
    t=t.drop(['pc','sib'],axis=0)

    return t


def cutoff_plotf(a,cnt,inb,hmmfold,outfold):

    models=["allLikelihoods.withdeg3_inphapProbs","allLikelihoods.nodeg3_inphapProbs","read_inppshap"]
    covs=[4,0.2,0.1,0.05]
    outf=outfold+"contam%s_inbred%s_model_performance_allroc_asc%s.csv.gz" %(cnt,inb,a)


    Tr=[]
    Fa=[]
    rl=[]
    cv=[]
    met=[]
    cut=[]

    for cov in covs:
        for model1 in models:
            for ct in np.linspace(0,3,100):


                if model1.split('_')[0]=="allLikelihoods.withdeg3":
                    model=model1.split('_')[0].split('.')[0] + '_inphapProbs'
                    un_flag='withdeg3'
                    table=hmmfold+"roc/contam%s/inbred%s/model_performance_%s/coverage%s/asc%s/filtered0_cut%s.csv.gz" %(cnt,inb,model,cov,a,ct)
                elif model1.split('_')[0]=="allLikelihoods.nodeg3":
                    model=model1.split('_')[0].split('.')[0] + '_inphapProbs'
                    un_flag='nodeg3'
                    table=hmmfold+"rocNo3deg/contam%s/inbred%s/model_performance_%s/coverage%s/asc%s/filtered0_cut%s.csv.gz" %(cnt,inb,model,cov,a,ct)
                elif model1.split('_')[0]=="read":
                    un_flag=0
                    model=model1
                    table="/mnt/diversity/divyaratan_popli/100arc/inbreeding/fastsim_gz_highDiversity/hmm_numba_corrA_relSim_fast_theoretical_sib_grA/rocNo3deg/contam%s/inbred%s/model_performance_%s/coverage%s/asc%s/filtered0_cut%s.csv.gz" %(cnt,inb,model,cov,a,ct)

                t= pd.read_csv(table, sep=",", header=0,index_col=0)
                t=all_manipulation(t)


                if un_flag=='nodeg3':
                    table_pcsib=hmmfold+"roc/contam%s/inbred%s/pc_sib_performance_%s/coverage%s/asc%s/filtered0_cut%s/pc_sib.csv.gz" %(cnt,inb,model,cov,a,ct)
                    tfir=pd.read_csv(table_pcsib, sep=",", header=0,index_col=0) #.loc[['pc','sib'],['pc','sib']]


                    print(tfir)

                    t=kin_manipulation(t)
                elif un_flag=='withdeg3':
                    t=kin_manipulation3(t)
                elif un_flag==0:
                    t=read_manipulation(t)

                if un_flag=='nodeg3' or un_flag==0:
                    for rel in ['un_all','sec','fir','id']:
                        TP=t.loc[rel,rel]/np.sum(t.loc[rel,:])
                        FP=(np.sum(t.loc[:,rel])-t.loc[rel,rel])/np.sum(t.loc[:,rel])
                        #plt.plot(ct,TP,marker='o',color='pink')
                        #plt.plot(ct,FP,marker='o',color='green')
                        Tr.append(TP)
                        Fa.append(FP)
                        rl.append(rel)
                        cv.append(cov)
                        met.append(model)
                        cut.append(ct)


                    if un_flag=='nodeg3':
                        for rel in ['pc','sib']:
                            TP=tfir.loc[rel,rel]/np.sum(tfir.loc[rel,:])
                            FP=(np.sum(tfir.loc[:,rel])-tfir.loc[rel,rel])/np.sum(tfir.loc[:,rel])

                            Tr.append(TP)
                            Fa.append(FP)
                            rl.append(rel)
                            cv.append(cov)
                            met.append(model)
                            cut.append(ct)

                elif un_flag=='withdeg3':
                    for rel in ['un_all','deg3']:
                        TP=t.loc[rel,rel]/np.sum(t.loc[rel,:])
                        FP=(np.sum(t.loc[:,rel])-t.loc[rel,rel])/np.sum(t.loc[:,rel])
                        if rel=='un_all':
                            rel='un_all_withdeg3'

                        Tr.append(TP)
                        Fa.append(FP)
                        rl.append(rel)
                        cv.append(cov)
                        met.append(model)
                        cut.append(ct)


    df=pd.DataFrame(
        {
         'True_positive': Tr,
         'False_positive': Fa,
         'Relatedness': rl,
         'Coverage': cv,
         'method': met,
         'cutoff': cut
         })

    df=df.dropna()
    with pd.option_context('display.max_rows', len(df.index), 'display.max_columns', len(df.columns)):
                df.to_csv(outf, sep=',')
    return df



def comparison_table(a,cnt,inb,hmmfold):

    models=["allLikelihoods.withdeg3_inphapProbs","allLikelihoods.nodeg3_inphapProbs","read_inppshap"]
    covs=[4,0.2,0.1,0.05]
    #outf=outfold+"contam%s_inbred%s_model_performance_allroc_asc%s.csv.gz" %(cnt,inb,a)


    Tr=[]
    Fa=[]
    rl=[]
    cv=[]
    met=[]
    cut=[]

    for cov in covs:
        for model1 in models:

            ct=1.0

            if model1.split('_')[0]=="allLikelihoods.withdeg3":
                model=model1.split('_')[0].split('.')[0] + '_inphapProbs'
                un_flag='withdeg3'
                table=hmmfold+"roc/contam%s/inbred%s/model_performance_%s/coverage%s/asc%s/filtered0_cut%s.csv.gz" %(cnt,inb,model,cov,a,ct)
            elif model1.split('_')[0]=="allLikelihoods.nodeg3":
                model=model1.split('_')[0].split('.')[0] + '_inphapProbs'
                un_flag='nodeg3'
                table=hmmfold+"rocNo3deg/contam%s/inbred%s/model_performance_%s/coverage%s/asc%s/filtered0_cut%s.csv.gz" %(cnt,inb,model,cov,a,ct)
            elif model1.split('_')[0]=="read":
                un_flag=0
                model=model1
                table="/mnt/diversity/divyaratan_popli/100arc/inbreeding/fastsim_gz_highDiversity/hmm_numba_corrA_relSim_fast_theoretical_sib_grA/rocNo3deg/contam%s/inbred%s/model_performance_%s/coverage%s/asc%s/filtered0_cut%s.csv.gz" %(cnt,inb,model,cov,a,ct)

            t= pd.read_csv(table, sep=",", header=0,index_col=0)
            t=all_manipulation(t)


            if un_flag=='nodeg3':
                table_pcsib=hmmfold+"roc/contam%s/inbred%s/pc_sib_performance_%s/coverage%s/asc%s/filtered0_cut%s/pc_sib.csv.gz" %(cnt,inb,model,cov,a,ct)
                tfir=pd.read_csv(table_pcsib, sep=",", header=0,index_col=0) #.loc[['pc','sib'],['pc','sib']]


                print(tfir)

                t=kin_manipulation(t)
            elif un_flag=='withdeg3':
                t=kin_manipulation3(t)
            elif un_flag==0:
                t=read_manipulation(t)

            if un_flag=='nodeg3' or un_flag==0:
                for rel in ['un_all','sec','fir','id']:
                    TP=t.loc[rel,rel]/np.sum(t.loc[rel,:])
                    FP=(np.sum(t.loc[:,rel])-t.loc[rel,rel])/np.sum(t.loc[:,rel])
                    #plt.plot(ct,TP,marker='o',color='pink')
                    #plt.plot(ct,FP,marker='o',color='green')
                    Tr.append(TP)
                    Fa.append(FP)
                    rl.append(rel)
                    cv.append(cov)
                    met.append(model)
                    cut.append(ct)


                if un_flag=='nodeg3':
                    for rel in ['pc','sib']:
                        TP=tfir.loc[rel,rel]/np.sum(tfir.loc[rel,:])
                        FP=(np.sum(tfir.loc[:,rel])-tfir.loc[rel,rel])/np.sum(tfir.loc[:,rel])

                        Tr.append(TP)
                        Fa.append(FP)
                        rl.append(rel)
                        cv.append(cov)
                        met.append(model)
                        cut.append(ct)

            elif un_flag=='withdeg3':
                for rel in ['un_all','deg3']:
                    TP=t.loc[rel,rel]/np.sum(t.loc[rel,:])
                    FP=(np.sum(t.loc[:,rel])-t.loc[rel,rel])/np.sum(t.loc[:,rel])
                    if rel=='un_all':
                        rel='un_all_withdeg3'

                    Tr.append(TP)
                    Fa.append(FP)
                    rl.append(rel)
                    cv.append(cov)
                    met.append(model)
                    cut.append(ct)


    df=pd.DataFrame(
        {
         'True_positive': Tr,
         'False_positive': Fa,
         'Relatedness': rl,
         'Coverage': cv,
         'method': met,
         'cutoff': cut
         })

    df=df.dropna()
    #with pd.option_context('display.max_rows', len(df.index), 'display.max_columns', len(df.columns)):
    #            df.to_csv(outf, sep=',')
    return df



def comparison_plotf(hmmfold,outfold, outs):
    dff=pd.DataFrame(
        columns=['True_positive', 'False_positive', 'Relatedness', 'Coverage', 'method',
           'cutoff', 'Ascertainment', 'Contamination', 'ROH']
        )

    for pl in range(4):
        if pl==0:
            a=0
            cnt=0
            inb=0
        elif pl==1:
            a=0
            cnt=1
            inb=0
        elif pl==2:
            a=0
            cnt=0
            inb=1
        elif pl==3:
            a=2
            cnt=0
            inb=0
        #lab1="KIn_contam=%s,asc=%s,roh=%s" %(cnt,a,inb)
        #lab2="READ_contam=%s,asc=%s,roh=%s" %(cnt,a,inb)
        df=comparison_table(a=a, cnt=cnt, inb=inb, hmmfold=hmmfold)
        df['Ascertainment']=a
        df['Contamination']=cnt
        df['ROH']=inb

        dff=dff.append(df)


    dff['contAscRoh']=dff['Contamination'].astype(str)+ '_' +dff['Ascertainment'].astype(str)+ '_' +dff['ROH'].astype(str)
    dff=dff.loc[dff['cutoff']==1.0,:]
    outff=outfold+'comparison_plot_data.csv.gz'
    with pd.option_context('display.max_rows', len(dff.index), 'display.max_columns', len(dff.columns)):
                    dff.to_csv(outff, sep=',')

    supf=dff
    supf.loc[supf['method']=='allLikelihoods_inphapProbs','method']='KIN'
    supf.loc[supf['method']=='read_inppshap','method']='READ'

    supf.loc[supf['Relatedness']=='un_all','Relatedness']='Unrelated w/o 3rd Degree'
    supf.loc[supf['Relatedness']=='un_all_withdeg3','Relatedness']='Unrelated'
    supf.loc[supf['Relatedness']=='deg3','Relatedness']='3rd Degree'
    supf.loc[supf['Relatedness']=='fir','Relatedness']='1st Degree'
    supf.loc[supf['Relatedness']=='sec','Relatedness']='2nd Degree'
    supf.loc[supf['Relatedness']=='id','Relatedness']='Identical'
    supf.loc[supf['Relatedness']=='pc','Relatedness']='Parent-Child'
    supf.loc[supf['Relatedness']=='sib','Relatedness']='Siblings'

    supf=supf.drop('contAscRoh', axis=1)
    supf=supf.drop('cutoff', axis=1)
    supf.columns=['True_positive', 'False_positive', 'Relatedness', 'Coverage', 'Method','Ascertainment', 'Contamination', 'ROH']

    with pd.option_context('display.max_rows', len(supf.index), 'display.max_columns', len(supf.columns)):
                    supf.to_csv(outs, sep=',',index=False)


def IBDstates(fold, outf, list_inds, runlist):

    accuracy=[]
    cov_all=[]
    run_all=[]
    pair_all=[]
    rel_all=[]


    for rel in list(list_inds.keys()):
        pair=list_inds[rel]
        for run in runlist:

            if rel=='gr':
                    rel1='grand'
            else:
                rel1=rel

            if rel == 'un':
                tr=np.ones(220)
            elif rel == 'pc':
                tr=np.ones(220)*0.75
            elif rel == 'id':
                tr=np.ones(220)*0.5
            elif not rel in ['un','pc','id']:
                trf=fold+'sanity_check/inbred0/run%s/coverage4/contam0/asc0/inputMode_hapProbs_fil0/%s_state_all.gz' %(run,rel1)
                tr=np.loadtxt(trf, dtype='float', delimiter = ",")

            for cov in [4,0.5,0.2,0.1,0.05,0.03]:
                vname=fold+"hmm_numba_viterbi_all_win/contam0/inbred0/run%s/coverage%s/filtered0/asc0/relation_%s_file/res_inphapProbs/pw_%s.csv.gz" %(run,cov,rel,pair)
                viterbi=np.loadtxt(vname, dtype='float', delimiter = ",")

                ac=sum(viterbi == tr)/len(viterbi)
                accuracy.append(ac)
                rel_all.append(rel)
                cov_all.append(cov)
                pair_all.append(pair)
                run_all.append(run)

    df=pd.DataFrame({
        'pair': pair_all,
        'rel': rel_all,
        'run': run_all,
        'accuracy': accuracy,
        'cov': cov_all
        })

    with pd.option_context('display.max_rows', len(df.index), 'display.max_columns', len(df.columns)):
                df.to_csv(outf, sep=',',index=False)
