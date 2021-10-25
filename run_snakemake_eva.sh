#!/bin/bash

#run with is_contam=0

snakemake all \
-w 60 \
--jobname "{rule}_{jobid}" \
-j 20 \
--cluster-config cluster.yaml \
--cluster "qsub -q all.q -l h_vmem={cluster.mem} -o autosnake/smk.out -e autosnake/smk.err -cwd -V -S /bin/bash" $* \
--config is_contam=0
