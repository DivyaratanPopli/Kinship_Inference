#!/bin/bash

#launch as:
#bash {this_script}.sh {path_to}/{LIST}
#where LIST=a list containing the absoluthe path to all the bam files that are
#to be used in the kinship analysis

echo launching preparation to KIN_snakemake
echo make sure to be in the dir where the snakemake pipeline will be launched
echo make sure that the original file names are formatted as:
echo {sample_name}{.yadayada}.bam
echo and that it is unique

x=0
PATHLIST=$1

mkdir -p bamfiles
mkdir -p bamfiles/stephane
rm targets.txt

sed 's/.*\///g' $PATHLIST > filenames.tmp
tot=$(cat filenames.tmp | wc -l)

echo creating the symlinks
for i in $(cat filenames.tmp | sed 's/\..*//g')
do
  export y=$(( ${x} + 1 ))
  export x=$y
  echo ${x}/${length} creating the symlink for ${i} 

  path=$(grep $i $PATHLIST)
  ln -s $path bamfiles/stephane/${i}.bam
  echo ${i} >> targets.txt
done

echo creating the other lists
cp targets.txt unrelated.txt
cp targets.txt identical.txt

cp /home/luca_traverso/scripts/allsites.bed .

echo done!
echo $(ls targets.txt) $(wc -l targets.txt)
echo $(ls unrelated.txt) $(wc -l unrelated.txt)
echo $(ls identical.txt) $(wc -l identical.txt)

ls -altrh bamfiles/stephane/*.bam > samples_analyzed.txt
echo a list of all the samples of this analysis is in samples_analyzed.txt
