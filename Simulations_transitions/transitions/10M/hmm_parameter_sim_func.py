

import numpy as np
import pandas as pd

def get_transition(sibstate, grstate, hsibstate, avustate, deg3state, deg4state, deg5state, sibAfile, grandAfile,
hsibAfile, avuAfile, deg3Afile, deg4Afile, deg5Afile, no_st):


    states=np.loadtxt(fname=sibstate, dtype='float', delimiter=',')
    stateg=np.loadtxt(fname=grstate, dtype='float', delimiter=',')
    stateh=np.loadtxt(fname=hsibstate, dtype='float', delimiter=',')
    statea=np.loadtxt(fname=avustate, dtype='float', delimiter=',')
    state3=np.loadtxt(fname=deg3state, dtype='float', delimiter=',')
    state4=np.loadtxt(fname=deg4state, dtype='float', delimiter=',')
    state5=np.loadtxt(fname=deg5state, dtype='float', delimiter=',')

    state_array=[states, stateg, stateh, statea, state3, state4, state5]
    out=[sibAfile,grandAfile,hsibAfile,avuAfile,deg3Afile,deg4Afile,deg5Afile]


    for rel in range(len(state_array)):
        state=state_array[rel]
        sibA=np.zeros((no_st,no_st))
        #print(len(state))
        for j in range(22):

            for k in range(int(len(state)/22)-1):
                i=(j*10)+k
                #print(i)
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


        df_A=pd.DataFrame(sibA)

        #print(df_A)

        with pd.option_context('display.max_rows', len(df_A.index), 'display.max_columns', len(df_A.columns)):
                    df_A.to_csv(out[rel], sep=',', header=False, index=False)



def avgA(siblistA, grandlistA, hsiblistA, avulistA, deg3listA, deg4listA, deg5listA, e1, e, l, sibAfinal,
    grandAfinal, hsibAfinal, avuAfinal, pcAfinal,idAfinal,unAfinal,deg3Afinal,deg4Afinal,deg5Afinal):

    Asib=np.zeros((l,l))
    Ag=np.zeros((l,l))
    Ahsib=np.zeros((l,l))
    Aavu=np.zeros((l,l))
    Adeg3=np.zeros((l,l))
    Adeg4=np.zeros((l,l))
    Adeg5=np.zeros((l,l))

    for i in range(len(siblistA)):
        sibmat=pd.read_csv(siblistA[i], sep=",", header=None,index_col=False)
        gmat=pd.read_csv(grandlistA[i], sep=",", header=None,index_col=False)
        hsibmat=pd.read_csv(hsiblistA[i], sep=",", header=None,index_col=False)
        avumat=pd.read_csv(avulistA[i], sep=",", header=None,index_col=False)
        deg3mat=pd.read_csv(deg3listA[i], sep=",", header=None,index_col=False)
        deg4mat=pd.read_csv(deg4listA[i], sep=",", header=None,index_col=False)
        deg5mat=pd.read_csv(deg5listA[i], sep=",", header=None,index_col=False)

        Asib=Asib+sibmat
        Ag=Ag+gmat
        Ahsib=Ahsib+hsibmat
        Aavu=Aavu+avumat
        Adeg3=Adeg3+deg3mat
        Adeg4=Adeg4+deg4mat
        Adeg5=Adeg5+deg5mat

    Asib=pd.DataFrame(Asib/(i+1))
    Ag=pd.DataFrame(Ag/(i+1))
    Ahsib=pd.DataFrame(Ahsib/(i+1))
    Aavu=pd.DataFrame(Aavu/(i+1))
    Adeg3=pd.DataFrame(Adeg3/(i+1))
    Adeg4=pd.DataFrame(Adeg4/(i+1))
    Adeg5=pd.DataFrame(Adeg5/(i+1))
    Apc=pd.DataFrame(np.array([[1-(2*e1),e1,e1],[1-(2*e1),e1,e1],[1-(2*e1),e1,e1]]))
    Aun=pd.DataFrame(np.array([[e1,1-(2*e1),e1],[e1,1-(2*e1),e1],[e1,1-(2*e1),e1]]))
    Aid=pd.DataFrame(np.array([[e1,e1,1-(2*e1)],[e1,e1,1-(2*e1)],[e1,e1,1-(2*e1)]]))

    modA=[Asib,Ag,Ahsib,Aavu,Adeg3,Adeg4,Adeg5]

    for r in range(len(modA)):
        modA[r]=modA[r]/modA[r].sum(axis=1)[:,None]

        if r!=0:
            modA[r].loc[2,:]=[2/3 - e/2, 1/3 - e/2, e]
            modA[r].loc[:2,2]=e
            modA[r].loc[0,:1]=modA[r].loc[0,:1]-e/2
            modA[r].loc[1,:1]=modA[r].loc[1,:1]-e/2

    [Asib,Ag,Ahsib,Aavu,Adeg3,Adeg4,Adeg5]=modA


    with pd.option_context('display.max_rows', len(Asib.index), 'display.max_columns', len(Asib.columns)):
                Asib.to_csv(sibAfinal, sep=',', header=False, index=False)

    with pd.option_context('display.max_rows', len(Ag.index), 'display.max_columns', len(Ag.columns)):
                Ag.to_csv(grandAfinal, sep=',', header=False, index=False)

    with pd.option_context('display.max_rows', len(Ahsib.index), 'display.max_columns', len(Ahsib.columns)):
                Ahsib.to_csv(hsibAfinal, sep=',', header=False, index=False)

    with pd.option_context('display.max_rows', len(Aavu.index), 'display.max_columns', len(Aavu.columns)):
                Aavu.to_csv(avuAfinal, sep=',', header=False, index=False)

    with pd.option_context('display.max_rows', len(Apc.index), 'display.max_columns', len(Apc.columns)):
                Apc.to_csv(pcAfinal, sep=',', header=False, index=False)

    with pd.option_context('display.max_rows', len(Aun.index), 'display.max_columns', len(Aun.columns)):
                Aun.to_csv(unAfinal, sep=',', header=False, index=False)

    with pd.option_context('display.max_rows', len(Aid.index), 'display.max_columns', len(Aid.columns)):
                Aid.to_csv(idAfinal, sep=',', header=False, index=False)

    with pd.option_context('display.max_rows', len(Adeg3.index), 'display.max_columns', len(Adeg3.columns)):
                Adeg3.to_csv(deg3Afinal, sep=',', header=False, index=False)

    with pd.option_context('display.max_rows', len(Adeg4.index), 'display.max_columns', len(Adeg4.columns)):
                Adeg4.to_csv(deg4Afinal, sep=',', header=False, index=False)

    with pd.option_context('display.max_rows', len(Adeg5.index), 'display.max_columns', len(Adeg5.columns)):
                Adeg5.to_csv(deg5Afinal, sep=',', header=False, index=False)
