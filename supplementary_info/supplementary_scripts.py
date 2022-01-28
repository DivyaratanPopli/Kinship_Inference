import pandas as pd

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
