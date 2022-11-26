# MSC Thesis

## Legacy Imaging Survey, DR9

### Configurations needed to use this code

- Install requirements found in `requirements.txt` within the appropriate folder
- Additionally, install the DESI library `desitarget`, per instructions at https://github.com/desihub/desitarget

## Simulations

### Configurations needed to use this code

- Download PINOCCHIO from https://github.com/pigimonaco/Pinocchio and follow compilation instructions from the INSTALLATION file. This will require the following modules to be loaded:
	- gsl/2.6
	- fftw/3.3.10
	- gcc/4.8.5
	- python/3.7.4 (or maybe even higher could work)
	- openmpi/4.1.4

- For the subhalo extraction code, need hdf5 via `module load hdf5/1.10.7`, which enables usage of `h5py`.

- For the galaxy sampling methods, need to install PyCosmo (`https://cosmo-gitlab.phys.ethz.ch/cosmo/PyCosmo`, branch `published`), UCat (`https://cosmo-gitlab.phys.ethz.ch/cosmo/ucat`) and UFig (`https://cosmo-gitlab.phys.ethz.ch/cosmo/ufig`). The `requirements.txt` covers these installations, but the repositories need to be cloned beforehand.
	- Note: there is a current issue with UFig, which requires `ivy-wfengine>=0.3.1` but no such version exists. To solve this, we modify the `setup.py` file within UFig to require `ivy-wfengine` instead of `ivy-wfengine>=0.3.1`.
