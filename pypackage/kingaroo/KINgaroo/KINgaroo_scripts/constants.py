#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 15:59:44 2022

@author: divyaratan_popli
"""

#CHRM=list(range(1,23))
contam_est='contam_est_pairwise.txt'
id_contam_est='id_contam_est.txt'
splitbams, bedfiles, hapProbs, hmm_param, hbdf, likf ='splitbams/','bedfiles/','hapProbs/','hmm_parameters/','hbd_results/','hbd_results/likelihoods/'
outdiff, outtotal, id_outdiff, id_outtotal = 'input_diffs_hmm.csv', 'input_total_hmm.csv', 'input_hbd_hmm_diffs.csv', 'input_hbd_hmm_total.csv'
