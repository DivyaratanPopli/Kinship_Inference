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
<p>Here optional inputs are shown in [].

-h: help<br>
-bam: Path to directory containing bamfiles with chromosomes represented by numbers 1,2,..,X,Y<br>
-bed: Path to tab-separated .bed file containing chromosome, reference and alternate alleles at all<br> &nbsp;&nbsp;&nbsp;&nbsp;available positions ([see example](example_files/bedfile.bed))<br>
-T: path to file ([see example file](example_files/targets.txt))containing list of all bamfiles to be used in the analysis<br>
-cnt: We provide three options for contamination correction:<br>
  &nbsp;&nbsp;&nbsp;&nbsp;0: No contamination correction<br>
  &nbsp;&nbsp;&nbsp;&nbsp;1: Contamination correction using divergence between the target population and contaminating population. We<br>
     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;implement a script that uses an indexed vcf.gz file with an individual each from<br>
     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;target and contaminating populations.<br>
  &nbsp;&nbsp;&nbsp;&nbsp;0<cnt<1: Contamination correction using divergence value given by cnt<br>
-c: Number of cores (by default: all available cores)<br>
-i: Size of genomic windows in int, Options:10000000, 1000000 (by default we use 10000000)<br>
-t: Minimum number of nonzero windows for a library to be included in estimation for p_0 (by default:10)<br>
-cest: File with contamination estimates (for contamination correction)<br>
-vcf.gz: Compressed and indexed vcf file for calculation of divergence between target and contaminating populations<br>
-tar: Name of individual from target population in vcf.gz<br>
-cont: Name of individual from contaminating population in vcf.gz<br>
-r: Enter 1 to estimate long ROH, 0 to skip (by default 1)<br>

# Running KIN
```
KIN [-h] -i  -o  -T  [-r] [-c] [-t] [-p]
```
-h help<br>
-i Path to the folder where you ran KINgaroo<br>
-o Output location<br>
-T Path to file containing list of all bamfiles to be used in the analysis<br>
-r Location of directory containing ROH estimates (by default: same as -i)<br>
-c Cores (by default: all available cores)<br>
-t Minimum number of sites in a window from which ROH estimates are reliable used (by default: 10)<br>
-p p_0 estimate given by user (by default: Estimated from the data)<br>

# Output
The final results are available in the file KIN_results.csv ([see example](example_files/KIN_results.csv))<br>
The output file has following columns:<br>
-Pair: Name of all pairs<br>
-Relatedness: Most likely relation<br>
-Second Guess: Second most likely relation (outside the degree of the most likely relation)<br>
-Log Likelihood Ratio: Log likelihood ratio for above mentioned relations<br>
-Within Degree Second Guess: Second most likely relation within the relatedness degree of most likely relation<br>
-Within Degree Log Likelihood Ratio:Log likelihood ratio for within-degree relations<br>
-k0: Proportion of genome with no IBD sharing<br>
-k1: Proportion of genome with one chromosome in IBD<br>
-k2: Proportion of genome with both chromosomes in IBD<br>
-IBD Length: Total number of windows in IBD<br>
-IBD Number: Total number of IBD segments<br>

We distinguish between the columns 'Second Guess' and 'Within Degree Second Guess' as well as between 'Log Likelihood Ratio'<br>
and 'Within Degree Log Likelihood Ratio'. This becomes important in case of classification to siblings or parent-child,<br> where we want to know how certain we are that the pair is first degree relative as indicated by 'Log Likelihood Ratio', but<br>
we also want to know the certainty associated with classification as parent-child compared to siblings or vice-versa.</p>
