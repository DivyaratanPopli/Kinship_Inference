for cv in 4 0.5 0.2 0.1 0.05 0.03
do
	for i in 0 1
	do
		for a in 0 2
		do
			for cnt in 0 1
			do
				for fold in allLikelihoods_inphapProbs
				do
					snakemake all --cores 55 --config inbs=$i cov=$cv asc=$a cnt=$cnt folder=$fold runs=$r --unlock
					snakemake all --cores 55 --config inbs=$i cov=$cv asc=$a cnt=$cnt folder=$fold runs=$r
					
				done
			done
		done
	done
done
