# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 10:44:57 2022

@author: divyaratan_popli
"""
import numpy as np

STATES=[0,1,2]
CHRM=list(range(1,23))
RELS=['un','deg5','deg4','deg3','gr','hsib','avu','sib','pc','id']

#transitions for 10M windows
A_un=np.array([[1.00000e-06, 9.99998e-01, 1.00000e-06],
       [1.00000e-06, 9.99998e-01, 1.00000e-06],
       [1.00000e-06, 9.99998e-01, 1.00000e-06]])

A_deg5=np.array([[6.25593827e-01, 3.74405173e-01, 1.00000000e-06],
       [2.51892991e-02, 9.74809701e-01, 1.00000000e-06],
       [6.66666167e-01, 3.33332833e-01, 1.00000000e-06]])

A_deg4=np.array([[7.53090860e-01, 2.46908140e-01, 1.00000000e-06],
       [3.58180318e-02, 9.64180968e-01, 1.00000000e-06],
       [6.66666167e-01, 3.33332833e-01, 1.00000000e-06]])

A_deg3=np.array([[8.27500554e-01, 1.72498446e-01, 1.00000000e-06],
       [5.86467863e-02, 9.41352214e-01, 1.00000000e-06],
       [6.66666167e-01, 3.33332833e-01, 1.00000000e-06]])

A_gr=np.array([[9.09364377e-01, 9.06346235e-02, 1.00000000e-06],
       [9.06346235e-02, 9.09364377e-01, 1.00000000e-06],
       [6.66665667e-01, 3.33333333e-01, 1.00000000e-06]])

A_hsib=np.array([[8.35731279e-01, 1.64267721e-01, 1.00000000e-06],
       [1.64387619e-01, 8.35611381e-01, 1.00000000e-06],
       [6.66666167e-01, 3.33332833e-01, 1.00000000e-06]])

A_avu=np.array([[8.05569361e-01, 1.94429639e-01, 1.00000000e-06],
       [1.94943269e-01, 8.05055731e-01, 1.00000000e-06],
       [6.66666167e-01, 3.33332833e-01, 1.00000000e-06]])

A_sib=np.array([[0.72466448, 0.13766776, 0.13766776],
       [0.27533552, 0.69749226, 0.02717222],
       [0.27533552, 0.02717222, 0.69749226]])

A_pc=np.array([[9.99998e-01, 1.00000e-06, 1.00000e-06],
       [9.99998e-01, 1.00000e-06, 1.00000e-06],
       [9.99998e-01, 1.00000e-06, 1.00000e-06]])

A_id=np.array([[1.00000e-06, 1.00000e-06, 9.99998e-01],
       [1.00000e-06, 1.00000e-06, 9.99998e-01],
       [1.00000e-06, 1.00000e-06, 9.99998e-01]])

A_ALL10=[A_un, A_deg5, A_deg4, A_deg3, A_gr, A_hsib, A_avu, A_sib, A_pc, A_id]

#transitions for 1M windows

A1_un=np.array([[1.00000e-06, 9.99998e-01, 1.00000e-06],
       [1.00000e-06, 9.99998e-01, 1.00000e-06],
       [1.00000e-06, 9.99998e-01, 1.00000e-06]])

A1_deg5=np.array([[9.44469251e-01, 5.55297488e-02, 1.00000000e-06],
       [3.79712288e-03, 9.96201877e-01, 1.00000000e-06],
       [6.66666167e-01, 3.33332833e-01, 1.00000000e-06]])

A1_deg4=np.array([[9.65172598e-01, 3.48264021e-02, 1.00000000e-06],
       [5.22537969e-03, 9.94773620e-01, 1.00000000e-06],
       [6.66666167e-01, 3.33332833e-01, 1.00000000e-06]])

A1_deg3=np.array([[9.75351446e-01, 2.46475544e-02, 1.00000000e-06],
       [8.89350472e-03, 9.91105495e-01, 1.00000000e-06],
       [6.66666167e-01, 3.33332833e-01, 1.00000000e-06]])

A1_gr=np.array([[9.90098337e-01, 9.90066335e-03, 1.00000000e-06],
       [9.90066335e-03, 9.90098337e-01, 1.00000000e-06],
       [6.66665667e-01, 3.33333333e-01, 1.00000000e-06]])

A1_hsib=np.array([[9.76757155e-01, 2.32418453e-02, 1.00000000e-06],
       [2.38486625e-02, 9.76150338e-01, 1.00000000e-06],
       [6.66666167e-01, 3.33332833e-01, 1.00000000e-06]])

A1_avu=np.array([[9.71820131e-01, 2.81788691e-02, 1.00000000e-06],
       [2.87937922e-02, 9.71205208e-01, 1.00000000e-06],
       [6.66666167e-01, 3.33332833e-01, 1.00000000e-06]])

A1_sib=np.array([[9.61558173e-01, 1.92209134e-02, 1.92209134e-02],
       [3.84418268e-02, 9.61173806e-01, 3.84367020e-04],
       [3.84418268e-02, 3.84367020e-04, 9.61173806e-01]])

A1_pc=np.array([[9.99998e-01, 1.00000e-06, 1.00000e-06],
       [9.99998e-01, 1.00000e-06, 1.00000e-06],
       [9.99998e-01, 1.00000e-06, 1.00000e-06]])

A1_id=np.array([[1.00000e-06, 1.00000e-06, 9.99998e-01],
       [1.00000e-06, 1.00000e-06, 9.99998e-01],
       [1.00000e-06, 1.00000e-06, 9.99998e-01]])

A_ALL1=[A1_un, A1_deg5, A1_deg4, A1_deg3, A1_gr, A1_hsib, A1_avu, A1_sib, A1_pc, A1_id]