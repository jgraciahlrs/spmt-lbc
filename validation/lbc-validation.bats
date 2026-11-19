#!/usr/bin/bash
setup () {
    bats_load_library 'bats-support'
    bats_load_library 'bats-assert'
    bats_load_library 'bats-file'
}

@test "lbc-bats present" {
    assert_file_exists "${BATS_TEST_DIRNAME}/lbc-bats"
}

@test "LBC 4 MPI processes" {
    MPIEXEC_OPTS="-n 4 -ppn 4 -depth 24 --cpu-bind depth --gpu-bind per_task:1 "
    cd "${BATS_TEST_TMPDIR}"
    cp ${BATS_TEST_DIRNAME}/4-mpi-lbc.params ./lbc.params
    assert_file_exists ./lbc.params
    run mpiexec ${MPIEXEC_OPTS} ${BATS_TEST_DIRNAME}/lbc-bats
    assert_success
    assert_output --partial "Total density:    8129015.981518 at timestep     10"

    tail -n 1 ${BATS_TEST_TMPDIR}/performance.res >> "${BATS_TEST_DIRNAME}/performance.res"
}

@test "LBC 1 MPI processes" {
    MPIEXEC_OPTS="-n 1 -ppn 4 -depth 24 --cpu-bind depth --gpu-bind per_task:1 "
    cd "${BATS_TEST_TMPDIR}"
    cp ${BATS_TEST_DIRNAME}/1-mpi-lbc.params ./lbc.params
    assert_file_exists ./lbc.params
    run mpiexec ${MPIEXEC_OPTS} ${BATS_TEST_DIRNAME}/lbc-bats
    assert_success
    assert_output --partial "Total density:    8129015.985189 at timestep     10"

    tail -n 1 ${BATS_TEST_TMPDIR}/performance.res >> "${BATS_TEST_DIRNAME}/performance.res"
}
