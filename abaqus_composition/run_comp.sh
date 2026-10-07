#!/usr/bin/env bash
# run_comp.sh -- Linux counterpart of run_comp.ps1 -Native (Keio server): build the
# composition UMAT with the point model in Fortran (make_umat.py --native, no material
# server) and run one Abaqus job.
#
#   abaqus_composition/run_comp.sh abaqus_composition/one_elem.inp [case] [cpus]
#
# - work dir $WORKROOT/comp_<job> (default $HOME/abaqus_work);
# - case constants from write_eco_cfg.py <case> into BIOFILM_ECO_CASE (default 2sp_case6);
# - mp_mode=threads: the UMAT/UMATHT/UEL share the nutrient field through a Fortran
#   module, which only works when all CPUs are threads of one process (not MPI ranks);
# - prints PASS / NOT COMPLETE from the .sta file, as run_comp.ps1 does.
# The server mode of run_comp.ps1 (C shim + material_server.py) is not ported.
set -euo pipefail
inp=${1:?usage: run_comp.sh <inp> [case] [cpus]}
case=${2:-2sp_case6}
cpus=${3:-1}
repo=$(cd "$(dirname "$0")/.." && pwd)
py=${PYTHON:-python3}
abq=${ABAQUS:-abaqus}
job=$(basename "$inp" .inp)
wd=${WORKROOT:-$HOME/abaqus_work}/comp_$job
rm -rf "$wd"; mkdir -p "$wd"
cp "$repo/$inp" "$wd/"
# make_cube_inp.py writes <job>.seed.txt next to the .inp while compare_ansys.py
# reads it next to the .dat: carry it over, or the comparison cannot find the seed.
seed="$repo/${inp%.inp}.seed.txt"
if [ -f "$seed" ]; then cp "$seed" "$wd/"; fi
"$py" "$repo/abaqus_composition/make_umat.py" "$wd/umat_comp.for" --native
"$py" "$repo/abaqus_composition/write_eco_cfg.py" "$case" "$wd/eco_case.txt"
export BIOFILM_ECO_CASE="$wd/eco_case.txt"
cd "$wd"
set +e
"$abq" job="$job" input="$job.inp" user=umat_comp.for cpus="$cpus" mp_mode=threads \
    interactive ask_delete=OFF 2>&1 | tee abaqus_out.txt | tail -n 15
set -e
if [ -f "$job.sta" ] && grep -q 'COMPLETED SUCCESSFULLY' "$job.sta"; then
    echo "PASS: $job completed (native point model)"
else
    echo "NOT COMPLETE: see $wd ($job.msg / $job.log)"
fi
echo "work dir: $wd"
