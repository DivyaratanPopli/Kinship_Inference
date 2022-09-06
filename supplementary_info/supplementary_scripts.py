import pandas as pd
import numpy as np



def suppTable2(resultsf,kinf,outf2,outf3):
    res = pd.read_excel(resultsf,usecols=(0,1,5,6,7,8,16,17),header=None)

    res.columns=list(res.loc[3,:])
    res=res.drop([0,1,2,3],axis=0)
    res.reset_index(drop=True, inplace=True)

    res['IND1']=res['IND1'].replace({'_':''}, regex=True)
    res['IND2']=res['IND2'].replace({'_':''}, regex=True)

    res['pair']=res['IND1'] + '_' + res['IND2']
    res.loc[0,'IND1']


    kin=pd.read_csv(kinf, sep="\t", header=0, index_col=0)
    kin.loc[0,:]

    kin.loc[kin['relatedness']=='sib','relatedness_r'] ='First Degree'
    kin.loc[kin['relatedness']=='pc','relatedness_r'] ='First Degree'
    kin.loc[kin['relatedness']=='sec','relatedness_r'] ='Second Degree'
    kin.loc[kin['relatedness']=='id','relatedness_r'] ='IdenticalTwins/SameIndividual'
    kin.loc[kin['relatedness']=='un','relatedness_r'] ='Unrelated'
    kin.loc[kin['relatedness']=='deg3','relatedness_r'] ='Third Degree'

    kin.loc[kin['relatedness']=='sib','relatedness_l'] ='Sibling'
    kin.loc[kin['relatedness']=='pc','relatedness_l'] ='Parent-offspring'
    kin.loc[kin['relatedness']=='sec','relatedness_l'] ='Second Degree'
    kin.loc[kin['relatedness']=='id','relatedness_l'] ='Identical'
    kin.loc[kin['relatedness']=='un','relatedness_l'] ='Unrelated'
    kin.loc[kin['relatedness']=='deg3','relatedness_l'] ='Third Degree'


    ares=pd.merge(res, kin, how='inner', on=['pair'])
    ares.loc[ares['Relationship (r/k0)* ']=='unrelated', 'Relationship (r/k0)* ']='Unrelated'
    ares.loc[ares['Relationship (r/k0)* ']=='2nd degree', 'Relationship (r/k0)* ']='Second Degree'

    #where do read and kin differ?

    diff=ares.loc[ares['relatedness_r']!=ares['Relationship (READ)'],['relatedness_r','relatedness_l','loglik_ratio','Relationship (READ)','Z_upper','Z_lower','Relationship (r/k0)* ','pair','nSNPs']]
    diff.reset_index(drop=True, inplace=True)


    dtable=diff.copy()
    dtable.columns=['relatedness_r', 'Relatedness KIn', 'Log liklihood ratio (KIn)', 'Relatedness READ','Z_upper (READ)', 'Z_lower (READ)', 'Relatedness lcMLkin', 'pair', 'nSNPs']
    dtable=dtable.drop('relatedness_r',1)

    #looking at cases where read says first degree and others disagree
    fread=ares.loc[ares['Relationship (READ)']=='First Degree',['relatedness_r','relatedness_l','withinDeg_ll','Relationship (READ)','Z_upper','Z_lower','Relationship (r/k0)* ','pair','nSNPs']]
    rtable=fread.copy()
    rtable.columns=['relatedness_r', 'Relatedness KIn', 'Log liklihood ratio (KIn)', 'Relatedness READ','Z_upper (READ)', 'Z_lower (READ)', 'Relatedness lcMLkin', 'pair', 'nSNPs']
    rtable=rtable.drop('relatedness_r',1)

    with pd.option_context('display.max_rows', len(dtable.index), 'display.max_columns', len(dtable.columns)):
                dtable.to_csv(outf2, sep=',',index=False)

    with pd.option_context('display.max_rows', len(rtable.index), 'display.max_columns', len(rtable.columns)):
                rtable.to_csv(outf3, sep=',',index=False)







def cnt_manipulation3(t1):

    t=t1.copy()
    t['un_all']=t['un']+t['deg5']+t['deg4']
    t=t.drop(['un','deg5','deg4'],axis=1)

    t['fir']=t['sib']+t['pc']
    t=t.drop(['pc','sib'],axis=1)

    t['sec']=t['gr']+t['avu']+t['hsib']
    t=t.drop(['gr','avu','hsib'],axis=1)


    t = t.reindex(index=['un_all', 'deg3', 'sec', 'fir', 'id'], columns=['un_all', 'deg3', 'sec', 'fir', 'id', 'NF'])

    t.loc['un_all',:]=t.loc['un_all',:] * 8
    t.loc['deg3',:]=t.loc['deg3',:] * 1
    t.loc['sec',:]=t.loc['sec',:] * 2
    t.loc['fir',:]=t.loc['fir',:] * 4
    t.loc['id',:]=t.loc['id',:] * 0

    return t

def cntall_manipulation(t1):
    t=t1.copy()
    t.loc['un_all',:]=t.loc['un',:]
    t=t.drop(['un','deg4','deg5'],axis=0)
    t.loc['sec',:]=(t.loc['gr',:] + t.loc['hsib',:]) / 2
    t=t.drop(['gr','hsib','avu'],axis=0)
    t.loc['fir',:]=(t.loc['pc',:]* 3 + t.loc['sib',:]* 1)/4
    t=t.drop(['pc','sib'],axis=0)

    return t


def misscontam_plotf(outf,covs,ctlist):



    Tr=[]
    Fa=[]
    rl=[]
    cv=[]
    met=[]
    cnt=[]

    for cov in covs:
            for ct in ctlist:

                table="/mnt/diversity/divyaratan_popli/review_sim/hmm_numba/miss_contam/output/missC%s_10/contam1/inbred0/model_performance_allLikelihoods_inphapProbs/coverage%s/asc0/filtered0_cut1.0.csv.gz" %(ct,cov)

                t= pd.read_csv(table, sep=",", header=0,index_col=0)
                t=cntall_manipulation(t)
                #for ind10 rem id
                t=t.drop(['id'],axis=0)


                table_pcsib="/mnt/diversity/divyaratan_popli/review_sim/hmm_numba/miss_contam/output/missC%s_10/contam1/inbred0/pc_sib_performance_allLikelihoods_inphapProbs/coverage%s/asc0/filtered0_cut1/pc_sib.csv.gz" %(ct,cov)
                tfir=pd.read_csv(table_pcsib, sep=",", header=0,index_col=0) #.loc[['pc','sib'],['pc','sib']]




                t=cnt_manipulation3(t)



                for rel in ['un_all','deg3','sec','fir']:
                    TP=t.loc[rel,rel]/np.sum(t.loc[rel,:])
                    FP=(np.sum(t.loc[:,rel])-t.loc[rel,rel])/np.sum(t.loc[:,rel])
                    #plt.plot(ct,TP,marker='o',color='pink')
                    #plt.plot(ct,FP,marker='o',color='green')
                    Tr.append(TP)
                    Fa.append(FP)
                    rl.append(rel)
                    cv.append(cov)
                    cnt.append(ct)



                for rel in ['pc','sib']:
                    TP=tfir.loc[rel,rel]/np.sum(tfir.loc[rel,:])
                    FP=(np.sum(tfir.loc[:,rel])-tfir.loc[rel,rel])/np.sum(tfir.loc[:,rel])

                    Tr.append(TP)
                    Fa.append(FP)
                    rl.append(rel)
                    cv.append(cov)
                    cnt.append(ct)


    df=pd.DataFrame(
        {
         'True_positive': Tr,
         'False_positive': Fa,
         'Relatedness': rl,
         'Coverage': cv,
         'contamination': cnt
         })

    df=df.dropna()
    with pd.option_context('display.max_rows', len(df.index), 'display.max_columns', len(df.columns)):
                df.to_csv(outf, sep=',')
    return df
