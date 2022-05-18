for cv in 4 0.5 0.2 0.1 0.05 0.03
do
	for a in 0 2
	do
		for i in 0 1
		do
			for cn in 0 1
			do

				#snakemake allt --cores 140 --config inbs=$i inb=$i cov=$cv asc=$a cnt=$cn folder=$fold runs=$r
				snakemake all --cores 140 --config inbs=$i inb=$i cov=$cv asc=$a cnt=$cn folder=$fold runs=$r
			done
		done
	done
done
