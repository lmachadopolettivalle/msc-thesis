for mass_cut in "2e12" "8e12"
do
	for quenching_time in "1.0" "2.0"
	do
		/cluster/home/lmachado/venv/bin/python submit_2D_sham.py --region BASS --num_z_bins 50 --num_mass_bins 60 --mass_cut $mass_cut --quenching_time $quenching_time --particle_count_pinocchio 2048
	done
done
