KIN is a Hidden-Markov-Model-based approach to identify identity-by-descent fragments and to
estimate degree of relatedness from ancient DNA data. KIN can accurately determine up to
3rd-degree relatives and differentiate between sibling and parent-child relationships with
as little as 0.05x coverage.


# Conda Environment
KIN requires Python 3.8+ and relies on a number of non-standard libraries. Here is
the list of these dependencies with the versions that we used:

- scipy (version 1.8.0)
- numpy (version 1.21.1)
- pandas (version 1.3.1)
- numba (version 0.55.1)


We recommend using a conda environment with all these dependencies:
```
conda create -n test1 python=3.8 scipy=1.8.0 numpy=1.21.1 pandas=1.3.1 numba=0.55.1
```
# Installation
After downloading or cloning pypackage from github, you can install KIN
by typing from the terminal which should install all necessary dependencies:
```
pip3 install _path_to_kin
```

# Running KIN
```
KIN [-h] -I  -O  -T  [-r] [-c] [-t] [-p] [-i]
```
-h: Help<br>
-I: Path to the folder where you ran KINgaroo<br>
-O: Output location for KIN<br>
-T: Path to file containing list of all bamfiles to be used in the analysis (should be same as that used in previous package)<br>
-r: Location of directory containing ROH estimates (by default: same as -I)<br>
-c: Cores (by default: all available cores)<br>
-t: Minimum number of sites in a window from which ROH estimates are reliable used (by default: 10)<br>
-p: p_0 estimate given by user (by default: Estimated from the data)<br>
-i: Size of genomic windows in int, use the same size as for KINgaroo (by default:10000000)<br>
# Output
The final results are available in the file KIN_results.csv ([see example file](example_files/KIN_results.csv))<br>
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
