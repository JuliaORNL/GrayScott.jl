# proceed with caution - the `cd` command is vulnerable to command injection
# usage: `./job_leconte.sh <dir_name>`
# setup:
# 1.) make a new directory for a new run and enter it `mkdir 003 && cd 003`
# 2.) copy and edit a settings file for the new run (be sure to edit the name of the output .bp file according to the changed variables) `cp ../001/settings-files.json . && vim settings-files.json`
# 3.) run the job `./job_leconte.sh 003`

(cd ${1:-.} && mpirun -n 1 julia --project=$GS_DIR $GS_DIR/gray-scott.jl settings-files.json | tee -a run.log)
