#! /bin/bash

###########
## Input ##
###########

ben_arch=SKL-SP

###################
## Configuration ##
###################


## Platform specific settings
if [[ "${ben_arch}" == "SKL-SP" ]]
then
    # declare -a metrics=("PAPI_VEC_SP" "PAPI_SP_OPS" "PAPI_TOT_CYC")
    # declare -a metrics=("PAPI_L1_TCM" "PAPI_LST_INS")
    # declare -a metrics=("PAPI_TOT_CYC" "PAPI_L2_TCM" "PAPI_L2_TCA")
    # declare -a metrics=("PAPI_TOT_CYC")
    declare -a metrics=("PAPI_LST_INS" "PAPI_TOT_CYC" "PAPI_L1_TCM" "PAPI_L2_TCA" "PAPI_L2_TCM" "PAPI_NATIVE_FP_ARITH:128B_PACKED_SINGLE" "PAPI_NATIVE_FP_ARITH:256B_PACKED_SINGLE" "PAPI_NATIVE_FP_ARITH:512B_PACKED_SINGLE")
    declare -a threads=(1 2 3 4 5 6 7 8 9 10)
else 
    echo ${ben_arch} "is not a valid architecture! Exiting..."
    exit
fi

icpc -g -O3 -qopenmp -ipo -mtune=skylake -xCORE-AVX512 -o gshf-mrqdt3 gshf-mrqdt3.cc marqfit.cc Event.cc

## Common file setup
out_dir=tau/scaling
exe="./gshf-mrqdt3"

export TAU_MAKEFILE=$TAU_DIR/x86_64/lib/Makefile.tau-icpc-papi-ompt-openmp
export TAU_CALLPATH_DEPTH=100
export TAU_UNWIND=1

tau_exe="tau_exec -ebs -T serial,icpc,papi,tbb "


for thr in "${threads[@]}"
do
    export OMP_NUM_THREADS=${thr}
    for metric in "${metrics[@]}"
    do    

        export TAU_METRICS=TIME,${metric}
        ${tau_exe} ${exe}
        mkdir ${out_dir}/${metric}_${thr}
        mv MULTI__* ${out_dir}/${metric}_${thr}
    done
done
    

