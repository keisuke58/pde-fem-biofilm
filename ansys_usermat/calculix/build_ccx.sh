#!/bin/sh
# Build CalculiX 2.21 with the Abaqus UMAT umat_biofilm_visco.f linked in
# (Debian/Ubuntu: apt install libspooles-dev libarpack2-dev liblapack-dev
# libblas-dev gfortran). Run from the repository root:
#   sh ansys_usermat/calculix/build_ccx.sh <dir with ccx_2.21.src.tar.bz2>
set -e
umat="$(pwd)/umat_biofilm_visco.f"
cd "$1"
tar xjf ccx_2.21.src.tar.bz2
cd CalculiX/ccx_2.21/src
cp "$umat" umat.f
printf '      implicit real*8(a-h,o-z)\n' > ABA_PARAM.INC
perl -0pi -e 's/^CFLAGS = .*$/CFLAGS = -w -O2 -I \/usr\/include\/spooles -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE -DNETWORKOUT/m;
  s/^FFLAGS = .*$/FFLAGS = -w -O2 -fallow-argument-mismatch/m;
  s/^LIBS = \\\n(?:.*\\\n)*.*\n/LIBS = -lspooles -larpack -llapack -lblas -lpthread -lm -lc\n/m;
  s/^(ccx_2.21: \$\(OCCXMAIN\) ccx_2.21.a) +\$\(LIBS\)/$1/m' Makefile
make -j4
