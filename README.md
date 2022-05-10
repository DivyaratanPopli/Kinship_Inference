KIN is a Hidden-Markov-Model-based approach to identify identity-by-descent fragments and to
estimate degree of relatedness from ancient DNA data. KIN can accurately determine up to
3rd-degree relatives and differentiate between sibling and parent-child relationships with
as little as 0.05x coverage.

KINgaroo is a software to generate input files for KIN from bamfiles. Optionally,
KINgaroo incorporates an adjustment for contamination, and an additional model to estimate the
location of long runs of homozygosity. This helps KIN to improve classification accuracy.

# Conda Environment
KIN and KINgaroo require Python 3.8+ and rely on a number of non-standard libraries. Here is
the list of these dependencies with the versions that we used:

- scipy (version 1.8.0)
- numpy (version 1.21.1)
- pandas (version 1.3.1)
- numba (version 0.55.1)
- pysam (version 0.19.0)
- pybedtools (version 0.9.0)

We recommend using a conda environment with all these dependencies:
```
conda create -n test1 python=3.8 scipy=1.8.0 numpy=1.21.1 pandas=1.3.1 numba=0.55.1 pysam=0.19.0 pybedtools=0.9.0
```
# Installation
After downloading or cloning pypackage from this repository, you can install KINgaroo
by typing from the terminal which should install all necessary dependencies:
```
pip3 install _path_to_kingaroo
```
Similarly, install kin:
```
pip3 install _path_to_kin
```

# Running KINgaroo
You can run KINgaroo from the terminal by typing:
```
 KINgaroo [-h] -bam  -bed  -T  -cnt  [-c] [-i] [-t] [-cest] [-vcf.gz] [-tar] [-cont] [-r]
```
Here optional inputs are shown in [].
```
-h help
-bam Path to directory containing bamfiles with chromosomes represented by numbers 1,2,..,X,Y
-bed Path to tab-separated .bed file containing chromosome, reference and alternate alleles at all 
     available positions [see example file](example_files/bedfile.bed)
-T path to file containing list of all bamfiles to be used in the analysis
- cnt We provide three options for contamination correction:
  0: No contamination correction
  1: Contamination correction using divergence between the target population and contaminating population. We
     implement a script that uses an indexed vcf.gz file with an individual each from target and contaminating
     populations.
  0<cnt<1 Contamination correction using divergence value given by cnt
-c Number of cores (by default: all available cores)
-i Size of genomic windows in int, Options:10000000, 1000000 (by default we use 10000000)
-t Minimum number of nonzero windows for a library to be included in estimation for p_0 (by default:10)
-cest File with contamination estimates (for contamination correction)
-vcf.gz Compressed and indexed vcf file for calculation of divergence between target and contaminating populations
-tar Name of individual from target population in vcf.gz
-cont Name of individual from contaminating population in vcf.gz
-r Enter 1 to estimate long ROH, 0 to skip (by default 1)
```
# Running KIN
```
-h help
-i Path to the folder where you ran KINgaroo
-o Output location
-T Path to file containing list of all bamfiles to be used in the analysis
-r Location of directory containing ROH estimates (by default: same as -i)
-c Cores (by default: all available cores)
-t Minimum number of sites in a window from which ROH estimates are reliable used (by default: 10)
-p p_0 estimate given by user (by default: Estimated from the data)
```

The final results are available in the file KIN_results.csv
