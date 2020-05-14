# FP_ARITH_INST_RETIRED.SCALAR_SINGLE
# FP_ARITH_INST_RETIRED.128B_PACKED_SINGLE
# FP_ARITH_INST_RETIRED.256B_PACKED_SINGLE
# measuring flops on leonhard cluster E5-2697 v4

# FP_ARITH_INST_RETIRED.SCALAR_DOUBLE
# FP_ARITH_INST_RETIRED.128B_PACKED_DOUBLE
# FP_ARITH_INST_RETIRED.256B_PACKED_DOUBLE
# INST_RETIRED.X87

# output of libpfm4/examples/checkevent
# Requested Event: FP_ARITH_INST_RETIRED.SCALAR_SINGLE
# Actual    Event: bdw_ep::FP_ARITH_INST_RETIRED:SCALAR_SINGLE:k=1:u=1:e=0:i=0:c=0:t=0:intx=0:intxcp=0
# PMU            : Intel Broadwell EP
# IDX            : 421527621
# Codes          : 0x5302c7
# Requested Event: FP_ARITH_INST_RETIRED.128B_PACKED_SINGLE
# Actual    Event: bdw_ep::FP_ARITH_INST_RETIRED:128B_PACKED_SINGLE:k=1:u=1:e=0:i=0:c=0:t=0:intx=0:intxcp=0
# PMU            : Intel Broadwell EP
# IDX            : 421527621
# Codes          : 0x5308c7
# Requested Event: FP_ARITH_INST_RETIRED.256B_PACKED_SINGLE
# Actual    Event: bdw_ep::FP_ARITH_INST_RETIRED:256B_PACKED_SINGLE:k=1:u=1:e=0:i=0:c=0:t=0:intx=0:intxcp=0
# PMU            : Intel Broadwell EP
# IDX            : 421527621
# Codes          : 0x5320c7
# Requested Event: FP_ARITH_INST_RETIRED.SCALAR_DOUBLE
# Actual    Event: bdw_ep::FP_ARITH_INST_RETIRED:SCALAR_DOUBLE:k=1:u=1:e=0:i=0:c=0:t=0:intx=0:intxcp=0
# PMU            : Intel Broadwell EP
# IDX            : 421527621
# Codes          : 0x5301c7
# Requested Event: FP_ARITH_INST_RETIRED.128B_PACKED_DOUBLE
# Actual    Event: bdw_ep::FP_ARITH_INST_RETIRED:128B_PACKED_DOUBLE:k=1:u=1:e=0:i=0:c=0:t=0:intx=0:intxcp=0
# PMU            : Intel Broadwell EP
# IDX            : 421527621
# Codes          : 0x5304c7
# Requested Event: FP_ARITH_INST_RETIRED.256B_PACKED_DOUBLE
# Actual    Event: bdw_ep::FP_ARITH_INST_RETIRED:256B_PACKED_DOUBLE:k=1:u=1:e=0:i=0:c=0:t=0:intx=0:intxcp=0
# PMU            : Intel Broadwell EP
# IDX            : 421527621
# Codes          : 0x5310c7
# Requested Event: FP_ARITH_INST_RETIRED.SCALAR
# Actual    Event: bdw_ep::FP_ARITH_INST_RETIRED:SCALAR_DOUBLE:SCALAR_SINGLE:k=1:u=1:e=0:i=0:c=0:t=0:intx=0:intxcp=0
# PMU            : Intel Broadwell EP
# IDX            : 421527621
# Codes          : 0x5303c7
# Requested Event: FP_ARITH_INST_RETIRED.PACKED
# Actual    Event: bdw_ep::FP_ARITH_INST_RETIRED:128B_PACKED_DOUBLE:128B_PACKED_SINGLE:256B_PACKED_DOUBLE:256B_PACKED_SINGLE:k=1:u=1:e=0:i=0:c=0:t=0:intx=0:intxcp=0
# PMU            : Intel Broadwell EP
# IDX            : 421527621
# Codes          : 0x533cc7
# Requested Event: FP_ARITH_INST_RETIRED.SINGLE
# Actual    Event: bdw_ep::FP_ARITH_INST_RETIRED:SCALAR_SINGLE:128B_PACKED_SINGLE:256B_PACKED_SINGLE:k=1:u=1:e=0:i=0:c=0:t=0:intx=0:intxcp=0
# PMU            : Intel Broadwell EP
# IDX            : 421527621
# Codes          : 0x532ac7
# Requested Event: FP_ARITH_INST_RETIRED.DOUBLE
# Actual    Event: bdw_ep::FP_ARITH_INST_RETIRED:SCALAR_DOUBLE:128B_PACKED_DOUBLE:256B_PACKED_DOUBLE:k=1:u=1:e=0:i=0:c=0:t=0:intx=0:intxcp=0
# PMU            : Intel Broadwell EP
# IDX            : 421527621
# Codes          : 0x5315c7

cmd="perf stat -e r5302c7,r5308c7,r5320c7,r5301c7,r5304c7,r5310c7,r5303c7,r533cc7,r532ac7,r5315c7
    $*"
echo $cmd
eval $cmd
