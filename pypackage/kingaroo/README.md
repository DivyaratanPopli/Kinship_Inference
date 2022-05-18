KINgaroo is a software to generate input files for KIN from processed bamfiles. Optionally,
KINgaroo incorporates an adjustment for contamination, and an additional model to estimate the
location of long runs of homozygosity. This helps KIN to improve classification accuracy.

# Conda Environment
KINgaroo requires Python 3.8+ and relies on a number of non-standard libraries. Here is
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

# Running KINgaroo
You can run KINgaroo from the terminal by typing:
```
  KINgaroo [-h] -bam  -bed  -T  -cnt  [-c] [-i] [-t] [-cest] [-d] [-tar] [-cont] [-r] [-p]
```
<p>Here optional inputs are shown in [].

-h: Help<br>
-bam: Path to directory containing bamfiles with chromosomes (represented by 1,2,..,X,Y) <br>
-bed: Path to tab-separated .bed file containing chromosome (1,2,..,X,..), reference and alternate alleles at all<br> &nbsp;&nbsp;&nbsp;&nbsp;available positions ([see example file](example_files/bedfile.bed))<br>
-T: Path to file ([see example file](example_files/targets.txt))containing list of all bamfiles to be used in the analysis<br>
-cnt: We provide three options for contamination correction:<br>
  &nbsp;&nbsp;&nbsp;&nbsp;0: No contamination correction<br>
  &nbsp;&nbsp;&nbsp;&nbsp;1: Contamination correction using divergence between the target population and contaminating population. Please<br>
     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;enter path to an indexed compressed vcf file [-d] with an individual each from<br>
     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;target [-tar] and contaminating populations [-cont]. Also required for this step: path<br>
     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;to contamination estimates file [-cest]<br>
  &nbsp;&nbsp;&nbsp;&nbsp;0<cnt<1: Contamination correction using divergence value entered here (0<cnt<1). Also required for this step: path<br>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;to contamination estimates file [-cest] <br>
-c: Number of cores (by default: all available cores)<br>
-i: Size of genomic windows in int, Options:10000000, 1000000 (by default we use 10000000)<br>
-t: Minimum number of nonzero windows for a library to be included in estimation for p_0 (by default:10)<br>
-cest: File with contamination estimates with 2 tab-separated columns: name,contamination<br>
-d: Compressed and indexed vcf file for calculation of divergence between target and contaminating populations. Please make sure that your vcf file<br>
 &nbsp;&nbsp;&nbsp;&nbsp;has genotypes [GT] represented in one of the following formats: X|Y (for phased files), X/Y (for unphased files),X (for pseudohaploids).<br>
 &nbsp;&nbsp;&nbsp;&nbsp;Here X,Y are 0/1 for ancestral/derived allele<br>
-tar: Name of individual from target population in [-d]<br>
-cont: Name of individual from contaminating population in [-d]<br>
-r: Enter 1 to estimate long ROH, 0 to skip (by default 1)<br>
-p: p_0 estimate given by user (by default: Estimated from the data)<br>
