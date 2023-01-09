for num_z_bins in "50" "75" "100" "125" "150" "175" "200" 
do
	for num_mass_bins in "30" "40" "50" "60" "100"
	do
		/cluster/home/lmachado/venv/bin/python submit_2D_sham.py --region BASS --num_z_bins $num_z_bins --num_mass_bins $num_mass_bins --mass_cut "8e12" --p_populated_subhalos 1 --p_red_subhalos 1 --particle_count_pinocchio 2048
	done
done
