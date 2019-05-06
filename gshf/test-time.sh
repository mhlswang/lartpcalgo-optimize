#! /bin/bash

###########
## Input ##
###########

ben_arch=SKL-SP

###################
## Configuration ##
###################

module load intel papi caliper

## Platform specific settings
if [[ "${ben_arch}" == "SKL-SP" ]]
then
    # declare -a metrics=("PAPI_TOT_CYC")
    declare -a threads=(1 2 3 4 5 6 7 8 9 10)
else 
    echo ${ben_arch} "is not a valid architecture! Exiting..."
    exit
fi

opt="-g -O3 -qopenmp -qopt-report=5 -mtune=skylake -xCORE-AVX512 "
icpc ${opt} -o gshf-mrqdt3 gshf-mrqdt3.cc marqfit.cc Event.cc

## Common file setup
exe="./gshf-mrqdt3"




for thr in "${threads[@]}"
do
    export OMP_NUM_THREADS=${thr}
    ${exe}
done
    

