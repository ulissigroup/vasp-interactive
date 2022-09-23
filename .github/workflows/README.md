# Running and maintaining the Github actions
------

## Unit test actions involved
- `package_and_unittest.yml`: Python lint and basic unit tests done on github action runner (VASP run with 1 core)
- `send_job_slurm.yaml`: Use scripts under the `tests/nersc_scriptes/` to enqueue a unit test job on corresponding slurm systems (Cori Haswell, Cori KNL, Perlmutter GPU etc).
- `*_status.yaml`: Action item only for showing the status of slurm unit test job (and update the README banner)

Note: 
- if ssh connection with NERSC is successful, `send_job_slurm.yaml` will always be passing
- At the end of the scripts under `tests/nersc_scriptes/`, the status signal will be sent to `*_status.yaml` and manually dispatch the action.

### Setting the slurm test environment
The steps to recreate the slurm test environment on NERSC Cori or Perlmutter is listed below. 
1. Create and activate a conda environment named `vpi` and install the basic dependencies `ase` `psutil` `pytest` 
2. Install the github cli by `conda install -c conda-forge gh`
3. Create a github token that has access to `vasp-interactive` repo (at least with read and actions privilege)
4. Login gh using credential `gh auth` and use the token created in step 3
5. Once successful, you should be able to see the `vasp-interactive` repo using `gh repo list`

Note: the above steps only need to be done once and can be used for both Cori and Perlmutter.



