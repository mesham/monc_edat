# MONC EDAT in-situ data analytics

This repository contains the MONC model with a port of the in-situ data analytics to use the EDAT task-based model, as well as the original in-situ data analytics code for comparison. 

## Pre-requisites
* NetCDF 4 
* FFTW
* FCM (https://www.metoffice.gov.uk/research/collaboration/fcm)
* MPICH
* EDAT

## Building MONC
MONC can either be built with GCC version 5 or the Cray compiler (8.4.1 tested.) To build you can use the FCM build system `fcm make -j4 -f fcm-make/monc-ubuntu-16.04-gnu-edat.cfg` for EDAT and `fcm make -j4 -f fcm-make/monc-ubuntu-16.04-gnu.cfg` for the original on a Linux system. For a Cray system then `fcm make -j4 -f fcm-make/monc-monc-cray-gnu-edat.cfg` and `fcm make -j4 -f fcm-make/monc-monc-cray-gnu.cfg` respectively. Note that these configuration files might need some modifications based on your exact system and set-up (see `fcm-make/env-ubuntu-16.04.cfg` and `env-cray.cfg` respectively.)

## Executing MONC
The built executable will be `build/bin/monc_driver_edat.exe` for the EDAT version and `build/bin/monc_driver.exe` for the original version. You can launch MONC with the `mpiexec` call, i.e. `mpiexec -np 2 ./build/bin/monc_driver_edat.exe --config=testcases/stratus/fire_sc.mcf` which will execute MONC on two processes using the stratus cloud test-case. 

There is also a PBS submission script for submission to larger scale machines

## Configuration
We have run this with the stratus test-case of testcases/stratus/fire_sc.mcf . Typically this will use `io/io_cfg_files/data_write_1file.xml` which is what we used for the benchmarking, however other smaller configurations are included for local runs (such as `io/io_cfg_files/edat_test.xml`) Additionally in the IO configuration directory there is a script, `duplicate_for_benchmarking.py` which will duplicate the data entries requested at a specific timestep, this enables us to request large amount of data without manually replicating lines in the file. Additionally there is an environment variable, `MONC_BENCHMARK_NUMS` which will duplicate all the data sent from the MONC computational core to the IO server for each request point any number of times (it is fine not to set this, the default is 1.)