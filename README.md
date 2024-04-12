# TacVar

## 1 Introduction


This repository is for validating the TacVar framework. Our developments have migrated to the PerfHound tool, constructing a more comprehensive visual system for hunting performance variability. Hence, this repo will only get documentation and bugfix updates..

TVkern generates run time variations which can be configured by given statistical distributions. By using different timing methods to measure the run time of the specific distribution, one can investigate the instability caused by timing fluctuations in parallel short-interval timing. Furthermore, based on the timing fluctuation distribution sample, TVfilt can help estimate the real distribution of timing results by decoupling the timing fluctuation from measurement results.


## 2 A simple use case

### 2.1 Preparing the test environment

### 2.2 Compiling source codes

### 2.3 Running VKern and visulizing timing fluctuations

### 2.4 Sampling timing fluctuations

### 2.5 Filtering the timing fluctuation with FilT

## 3 Options

## 3.1 TVkern

Using predefined macros (-D) to configure TVkern in compiling.

Mandatory options:
- USE_*: Available timing methods:
  - TSC: Serialized RDTSC/RDTSCP for x86 processors.
  - CGT: Call clock_gettime() using Linux monotonic clock.
  - PAPI: Call PAPI_get_real_nsec() in PAPI library.
  - PAPIX: In addition to PAPI, read 6 more event counters.
  - LIKWID: Call LIKWID marker APIS.
- V1:
- V2:
- FSIZE: Cache flush size in Bytes.

Optional options:
- NWARM: Time (ms) for warm up. (Default: 1000).
- NPASS: Ignoring the first NPASS test results. (Default: 2).
- NTEST: 
- TBASE:
- LCUT:
- HCUT:
- NPRECALC:

## 3.2 TVfilt

``` bash

$ ./filt.x --help
Usage: filt.x [OPTION...]

  -g, --hz-ns=FREQ           Tick per ns of the clock
  -l, --plow=PROB            Lowest threshold of possibility of a data bin
  -m, --met-file=FILE        Input measurement file
  -n, --nsamp=NSAMP          Number of samples in each optimization step
  -s, --sample-file=FILE     Input timing fluctuation file
  -w, --width=TIME           The least interval of a time bin (ns).
  -x, --cut-x=PROB           Cut the highest probability of met array
  -y, --cut-y=PROB           Cut the highest probability of timing fluctuation
                             array
  -?, --help                 Give this help list
      --usage                Give a short usage message

```