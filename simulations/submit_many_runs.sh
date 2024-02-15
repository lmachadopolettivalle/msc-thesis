for mass_cut in "6e12"
do
	for quenching_time in "2.0"
	do
		for region in "fullsky"
		do
			for seed in "666666"
			do
				/cluster/home/lmachado/venv/bin/python submit_2D_sham.py --region $region --num_z_bins 50 --num_mass_bins 60 --mass_cut $mass_cut --quenching_time $quenching_time --particle_count_pinocchio 2048 --seed $seed
			done
		done
	done
done
