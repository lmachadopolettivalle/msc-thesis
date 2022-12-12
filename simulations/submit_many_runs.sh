for mass_cut in "4e12" "6e12" "8e12" "1e13" "2e13" "4e13"
do
	/cluster/home/lmachado/venv/bin/python submit_2D_sham.py --region BASS --num_z_bins 150 --num_mass_bins 30 --mass_cut $mass_cut --p_populated_subhalos 1 --p_red_subhalos 1 --particle_count_pinocchio 512
done
