import numpy as np
import pandas as pd

FILENAME = "/cluster/home/lmachado/msc-thesis/simulations/explored_parameter_space.txt"

EXPLORED = pd.read_csv(FILENAME, sep=',', lineterminator='\n', index_col=None)

DEFAULT_NUM_Z_BINS = 150
DEFAULT_NUM_MASS_BINS = 30
DEFAULT_MASS_CUT = 8.0e12
DEFAULT_P_POPULATED_SUBHALOS = 1.0
DEFAULT_P_RED_SUBHALOS = 1.0

def get_new_ID():
    # Find an ID that does not yet exist in the DataFrame
    return 1 + max(EXPLORED["ID"])

def get_id_of_run(num_z_bins=DEFAULT_NUM_Z_BINS, num_mass_bins=DEFAULT_NUM_MASS_BINS, mass_cut=DEFAULT_MASS_CUT, p_populated_subhalos=DEFAULT_P_POPULATED_SUBHALOS, p_red_subhalos=DEFAULT_P_RED_SUBHALOS):
    run_ids = EXPLORED.loc[
        (EXPLORED["num_z_bins"] == num_z_bins) &
        (EXPLORED["num_mass_bins"] == num_mass_bins) &
        (np.isclose(EXPLORED["mass_cut"], mass_cut)) &
        (np.isclose(EXPLORED["p_populated_subhalos"], p_populated_subhalos)) &
        (np.isclose(EXPLORED["p_red_subhalos"], p_red_subhalos))
    ]["ID"].values

    if len(run_ids) == 0:
        return None

    return run_ids[0]

def store_new_run(num_z_bins=DEFAULT_NUM_Z_BINS, num_mass_bins=DEFAULT_NUM_MASS_BINS, mass_cut=DEFAULT_MASS_CUT, p_populated_subhalos=DEFAULT_P_POPULATED_SUBHALOS, p_red_subhalos=DEFAULT_P_RED_SUBHALOS):
    global EXPLORED

    # Make sure this run does not already exist
    ID = get_id_of_run(
        num_z_bins=num_z_bins,
        num_mass_bins=num_mass_bins,
        mass_cut=mass_cut,
        p_populated_subhalos=p_populated_subhalos,
        p_red_subhalos=p_red_subhalos,
    )
    if ID is not None:
        print(f"Run already exists with ID {ID}. Will not store a new entry.")
        return

    # Add run entry to DataFrame, and overwrite txt file
    new_ID = get_new_ID()
    new_df = pd.DataFrame(
        {
            "ID": [new_ID],
            "num_z_bins": [num_z_bins],
            "num_mass_bins": [num_mass_bins],
            "mass_cut": [mass_cut],
            "p_populated_subhalos": [p_populated_subhalos],
            "p_red_subhalos": [p_red_subhalos],
        }
    )

    print(f"Will store new run with ID {new_ID}.")

    EXPLORED = pd.concat(
        [EXPLORED, new_df],
        ignore_index=True,
    )

    EXPLORED.reset_index()

    EXPLORED.to_csv(FILENAME, index=False)

    return

def get_details_of_run(run_id):
    row = EXPLORED.loc[
        (EXPLORED["ID"] == run_id)
    ]
    if len(row) == 0:
        raise ValueError(f"Requested run ID {run_id} which is not present in the database.")

    data = row.to_dict(orient="list")

    data = {
        k: v[0]
        for k, v in data.items()
    }

    return data
