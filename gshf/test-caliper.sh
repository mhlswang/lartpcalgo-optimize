#! /bin/bash

###########
## Input ##
###########

ben_arch=SKL-SP

###################
## Configuration ##
###################

module load intel papi caliper

export CALI_CONFIG_FILE=~/soft/src/Caliper/examples/configs/papi_cycles.conf

## Platform specific settings
if [[ "${ben_arch}" == "SKL-SP" ]]
then
    # declare -a metrics=("PAPI_TOT_CYC")
    declare -a metrics=("PAPI_SP_OPS" "PAPI_TOT_CYC" "FP_ARITH:128B_PACKED_SINGLE" "FP_ARITH:256B_PACKED_SINGLE" "FP_ARITH:512B_PACKED_SINGLE" "FP_ARITH:SCALAR_SINGLE")
    # declare -a metrics=("PAPI_LST_INS" "PAPI_RES_STL" "PAPI_MEM_WCY" "PAPI_TOT_INS" "PAPI_SP_OPS" "PAPI_TOT_CYC" "L2_RQSTS:REFERENCES" "PAPI_L2_TCM" "PAPI_L1_TCM" "FP_ARITH:128B_PACKED_SINGLE" "FP_ARITH:256B_PACKED_SINGLE" "FP_ARITH:512B_PACKED_SINGLE" "FP_ARITH:SCALAR_SINGLE")
    declare -a threads=(1)
    # declare -a threads=(1 2 4 6 8 10)
else 
    echo ${ben_arch} "is not a valid architecture! Exiting..."
    exit
fi

opt="-g -O3 -qopenmp -qopt-report=5 -xSSE2 -DUSE_CALI -I${CALIPER_DIR}/include"
# opt="-g -O3 -qopenmp -qopt-report=5 -mtune=skylake -qopt-zmm-usage=high -xSKYLAKE-AVX512 -DUSE_CALI -I${CALIPER_DIR}/include"
lib="-L${CALIPER_DIR}/lib64 -lcaliper"
icpc ${opt} -o gshf-mrqdt3 gshf-mrqdt3.cc marqfit.cc Event.cc ${lib}

## Common file setup
out_dir=cali_sse_pragma
mkdir ${out_dir}
exe="./gshf-mrqdt3"




for thr in "${threads[@]}"
do
    export OMP_NUM_THREADS=${thr}
    for metric in "${metrics[@]}"
    do    

        export CALI_REPORT_FILENAME=${out_dir}/cali_${metric}_NTH${thr}.json
        export CALI_PAPI_COUNTERS=${metric}
        ${exe}

    done
done
    

