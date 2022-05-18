for cnt in 0 1
do
	for i in 0 1
	do
		for a in 0 2 3
		do
			for inp in hapProbs pshap
			do
				echo "$i"
				snakemake alltables --config i=$i a=$a inp=$inp cnt=$cnt --cores 120
			done
		done
	done
done
