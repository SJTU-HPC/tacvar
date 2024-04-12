#define _GNU_SOURCE
#define _ISOC11_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <stdlib.h>
#include <unistd.h>
#include <sched.h>
#include "mpi.h"

#ifdef TIMING
#if defined(USE_PAPI) || defined(USE_PAPIX6)
#include "papi.h"

#elif defined(USE_LIKWID)
#include "likwid-marker.h"
#define NEV 20
#endif

#endif

// Warmup for 1000ms.
#ifndef NWARM
#define NWARM 1000
#endif

// Number of tests for each interval
#ifndef NTEST
#define NTEST 10
#endif

#ifndef NARR
#define NARR 100
#endif

#ifndef NPASS
#define NPASS 1 
#endif


// Timing macros
#define _read_ns(_ns) \
    do {                                                \
        register uint64_t ns;                           \
        asm volatile(                                   \
            "\n\tRDTSCP"                                \
            "\n\tshl $32, %%rdx"                        \
            "\n\tor  %%rax, %%rdx"                      \
            "\n\tmov %%rdx, %0"                         \
            "\n\t"                                      \
            :"=r" (ns)                                  \
            :                                           \
            : "memory", "%rax", "%rdx");                \
        _ns = ns;                                       \
    } while(0);

#define _read_cy(_cy) \
    do {                                                \
        register uint64_t cy;                           \
        asm volatile(                                   \
            _pfc_read(0x40000001)                       \
            : "=r"(cy)                                  \
            :                                           \
            : "memory", "rax", "rcx", "rdx"             \
        );                                              \
        _cy = cy;                                       \
    } while(0)

#define _mfence asm volatile("lfence"   "\n\t":::);


#define NS_PER_TICK  1

void
tsc_start(uint64_t *cycle) {
    unsigned ch, cl;

    asm volatile (  "CPUID" "\n\t"
                    "RDTSC" "\n\t"
                    "mov %%edx, %0" "\n\t"
                    "mov %%eax, %1" "\n\t"
                    : "=r" (ch), "=r" (cl)
                    :
                    : "%rax", "%rbx", "%rcx", "%rdx");
    *cycle = ( ((uint64_t)ch << 32) | cl );
}

void
tsc_stop(uint64_t *cycle) {
    unsigned ch, cl;

    asm volatile (  "RDTSCP" "\n\t"
                    "mov %%edx, %0" "\n\t"
                    "mov %%eax, %1" "\n\t"
                    "CPUID" "\n\t"
                    : "=r" (ch), "=r" (cl)
                    :
                    : "%rax", "%rbx", "%rcx", "%rdx");

    *cycle = ( ((uint64_t)ch << 32) | cl );
}


/**
 * @brief Fill arr[size] with random number.
 */
void
fill_random(double *arr, size_t size) {
    struct timespec tv;
    uint64_t sec, nsec;
    clock_gettime(CLOCK_MONOTONIC, &tv);
    sec = tv.tv_sec;
    nsec = tv.tv_nsec;
    nsec = sec * 1e9 + nsec + NWARM * 1e6;
    srand(nsec);
    for (size_t i = 0; i < size; i ++){
        arr[i] = (float)rand() / (float)RAND_MAX;
    }

    return;
}

int
main(int argc, char **argv) {
    uint64_t ntest;
    /* Vars for Stencil */
    double **x, **y;
    double a = 0.21, b = 0.20;
    uint64_t narr;
    struct timespec tv;
    uint64_t volatile nsec_st, nsec_en; // For warmup
    int myrank, nrank, errid;
    double tsc_ns;

    if (argc >= 2) {
        narr = (uint64_t)atoll(argv[1]);
    } else {
        narr = NARR;
    }

    if (argc >= 3) {
        tsc_ns = atof(argv[2]);
    } else {
        tsc_ns = 1.0;
    }

    x = (double **)malloc(narr * sizeof(double*));
    y = (double **)malloc(narr * sizeof(double*));
    for (size_t i = 0; i < narr; i ++) {
        x[i] = (double *)malloc(narr * sizeof(double));
    }
    for (size_t i = 0; i < narr; i ++) {
        y[i] = (double *)malloc(narr * sizeof(double));
    }

    ntest = NPASS + NTEST;

    errid = MPI_Init(NULL, NULL);
    if (errid != MPI_SUCCESS) {
        printf("Faild to init MPI.\n");
        exit(1);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    
    if (myrank == 0) {
        printf("A 2D 5-point Jacobi stencil scheme.\nNTEST=%lu, NPASS=%u, NARR=%lu \n", 
                ntest, NPASS, narr);
    }

    // Warm up
    do {
        for (uint64_t i = 0; i < narr; i ++) {
            fill_random(x[i], narr);
            fill_random(y[i], narr);
        }
        if (myrank == 0) {
            printf("Warming up for %d ms.\n", NWARM);
        }
        struct timespec tv;
        uint64_t volatile sec, nsec; // For warmup
        clock_gettime(CLOCK_MONOTONIC, &tv);
        sec = tv.tv_sec;
        nsec = tv.tv_nsec;
        nsec = sec * 1e9 + nsec + NWARM * 1e6;

        while (tv.tv_sec * 1e9 + tv.tv_nsec < nsec) {
            clock_gettime(CLOCK_MONOTONIC, &tv);
            for (uint64_t i = 1; i < narr-1; i ++) {
                for (uint64_t j = 1; j < narr-1; j ++) {
                    y[i][j] = a * x[i][j] + b * (x[i-1][j] + x[i+1][j] + x[i][j-1] + x[i][j+1]);
                }
            }
        }
    } while (0);



#ifdef TIMING
    uint64_t *p_ns, ns0 = 0, ns1 = 0;

#ifdef USE_PAPI
    // Init PAPI
    int eventset = PAPI_NULL;
    PAPI_library_init(PAPI_VER_CURRENT);
    PAPI_create_eventset(&eventset);
    PAPI_start(eventset);
    
#elif USE_PAPIX6
    // Init PAPI
    int eventset = PAPI_NULL;
    int nev = 6;
    long long int ev_vals_0[6]={0}, ev_vals_1[6]={0};
    int64_t *p_ev;
    p_ev = (int64_t *)malloc(ntest * narr * nev * sizeof(int64_t));
    for (int iev = 0; iev < nev; iev ++) {
        ev_vals_0[iev] = 0;
        ev_vals_1[iev] = 0;
    }
    PAPI_library_init(PAPI_VER_CURRENT);
    PAPI_create_eventset(&eventset);
    PAPI_add_named_event(eventset, "cpu-cycles");
    PAPI_add_named_event(eventset, "instructions");
    PAPI_add_named_event(eventset, "cache-references");
    PAPI_add_named_event(eventset, "cache-misses");
    PAPI_add_named_event(eventset, "branches");
    PAPI_add_named_event(eventset, "branch-misses");
    PAPI_start(eventset);

#elif USE_LIKWID
    // Init LIKWID
    double ev_vals_0[NEV]={0}, ev_vals_1[NEV]={0}, time=0;
    int nev = NEV, count;
    int64_t *p_ev;
    LIKWID_MARKER_INIT;
    LIKWID_MARKER_THREADINIT;
    LIKWID_MARKER_REGISTER("vkern"); 
    LIKWID_MARKER_REGISTER("nev_count"); 
    LIKWID_MARKER_START("nev_count");
    LIKWID_MARKER_STOP("nev_count");
    LIKWID_MARKER_GET("nev_count", &nev, (double *)ev_vals_1, &time, &count);
    if (myrank == 0) {
        printf("LIKWID event count = %d\n", nev);
    }
    p_ev = (int64_t *)malloc(ntest * narr * nev * sizeof(int64_t));
    for (int iev = 0; iev < nev; iev ++) {
        ev_vals_0[iev] = 0;
        ev_vals_1[iev] = 0;
    }

#endif

    p_ns = (uint64_t *)malloc(ntest * narr * sizeof(uint64_t));
#endif

    for (uint64_t i = 0; i < narr; i ++) {
        fill_random(x[i], narr);
        fill_random(y[i], narr);
    }

    if (myrank == 0) {
        printf("Start running.\n");
        fflush(stdout);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    clock_gettime(CLOCK_MONOTONIC, &tv);
    nsec_st = tv.tv_sec * 1e9 + tv.tv_nsec;
    for (int it = 0; it < ntest; it ++) {
        for (uint64_t i = 0; i < narr; i ++) {
            for (uint64_t j = 0; j < narr; j ++) {
                x[i][j] = y[i][j];
            }
        }

        for (uint64_t j = 1; j < narr-1; j ++) {
#ifdef TIMING

// Timing.
#ifdef USE_PAPI
            ns0 = PAPI_get_real_nsec();

#elif USE_PAPIX6
            ns0 = PAPI_get_real_nsec();
            PAPI_read(eventset, ev_vals_0);

#elif USE_CGT
            clock_gettime(CLOCK_MONOTONIC, &tv);
            ns0 = tv.tv_sec * 1e9 + tv.tv_nsec;

#elif USE_WTIME
            ns0 = (uint64_t)(MPI_Wtime() * 1e9);

#elif USE_LIKWID
            //ns0 = 0;
            LIKWID_MARKER_START("vkern"); 

#elif USE_TSC
            tsc_start(&ns0);

#else
            _read_ns (ns0);
            _mfence;

#endif

#endif
            for (uint64_t k = 1; k < narr-1; k ++) {
                y[j][k] = a * x[j][k] + b * (x[j-1][k] + x[j+1][k] + x[j][k-1] + x[j][k+1]);
            }
#ifdef TIMING

#ifdef USE_PAPI
            ns1 = PAPI_get_real_nsec();
            p_ns[it*narr+j] = (uint64_t)(ns1 - ns0);

#elif USE_PAPIX6
            ns1 = PAPI_get_real_nsec();
            PAPI_read(eventset, ev_vals_1);
            for (int iev = 0; iev < nev; iev ++) {
                p_ev[it * narr * nev + j * nev + iev] = (int64_t)(ev_vals_1[iev] - ev_vals_0[iev]);
            }
            p_ns[it*narr+j] = (uint64_t)(ns1 - ns0);

#elif USE_CGT
            clock_gettime(CLOCK_MONOTONIC, &tv);
            ns1 = tv.tv_sec * 1e9 + tv.tv_nsec;
            p_ns[it*narr+j] = ns1 - ns0;

#elif USE_WTIME
            ns1 = (uint64_t)(MPI_Wtime() * 1e9);
            p_ns[it*narr+j] = ns1 - ns0;

#elif USE_LIKWID
            LIKWID_MARKER_STOP("vkern"); 
            LIKWID_MARKER_GET("vkern", &nev, (double*)ev_vals_1, &time, &count);
            for (int iev = 0; iev < nev; iev ++) {
                p_ev[it * narr * nev + j * nev + iev] = (int64_t)ev_vals_1[iev] - (int64_t)ev_vals_0[iev];
                ev_vals_0[iev] = ev_vals_1[iev];
            }
            // We do not use "time" argument as the timestamp because the perfmon swith the timer
            // unexpectedly. It is good to use FIXC2: CPU_CLK_UNHALTED_REF and convert with tsc_ns
            // ns1 = (uint64_t)((double)p_ev[it * narr * nev + j * nev + 2] / tsc_ns);
            ns1 = (uint64_t) (time * 1e9);
            p_ns[it*narr+j] = ns1 - ns0;
            ns0 = ns1;

#elif USE_TSC
            tsc_stop(&ns1);
            p_ns[it*narr+j] = ns1 - ns0;
            p_ns[it*narr+j] = (uint64_t)((double)(ns1 - ns0) / tsc_ns);

#else
            _read_ns (ns1);
            _mfence;
            p_ns[it*narr+j] = (uint64_t)((double)(ns1 - ns0) / tsc_ns);
#endif

#endif
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &tv);
    nsec_en = tv.tv_sec * 1e9 + tv.tv_nsec;
    printf("Rank %d run time: %lu ns\n", myrank, nsec_en - nsec_st);

    MPI_Barrier(MPI_COMM_WORLD);

#if defined(USE_PAPI) || defined(USE_PAPIX6)
    PAPI_shutdown();

#elif USE_LIKWID
    // Finalize LIKWID
    LIKWID_MARKER_CLOSE;    
#endif

    // Each rank writes its own file.
#ifdef TIMING
    char fname[4096], myhost[1024];
    gethostname(myhost, 1024);
#ifdef USE_PAPI
    sprintf(fname, "jacobi2d5p_papi_time_%d_%s.csv", myrank, myhost);
    FILE *fp = fopen(fname, "w");
#elif USE_CGT
    sprintf(fname, "jacobi2d5p_cgt_time_%d_%s.csv", myrank, myhost);
    FILE *fp = fopen(fname, "w");
#elif USE_WTIME
    sprintf(fname, "jacobi2d5p_wtime_time_%d_%s.csv", myrank, myhost);
    FILE *fp = fopen(fname, "w");
#elif USE_PAPIX6
    sprintf(fname, "jacobi2d5p_papix6_time_%d_%s.csv", myrank, myhost);
    FILE *fp = fopen(fname, "w");
#elif USE_LIKWID
    sprintf(fname, "jacobi2d5p_likwid_time_%d_%s.csv", myrank, myhost);
    FILE *fp = fopen(fname, "w");
#elif USE_TSC
    sprintf(fname, "jacobi2d5p_tsc_time_%d_%s.csv", myrank, myhost);
    FILE *fp = fopen(fname, "w");
#else
    sprintf(fname, "jacobi2d5p_stiming_time_%d_%s.csv", myrank, myhost);
    FILE *fp = fopen(fname, "w");
#endif

    for (int it = NPASS; it < ntest; it ++) {
        for (size_t j = 1; j < narr-1; j ++) {
            fprintf(fp, "%d,%lu", myrank, p_ns[it*narr+j]);

#if defined(USE_LIKWID) || defined(USE_PAPIX6)
            for (int iev = 0; iev < nev; iev ++) {
                fprintf(fp, ",%ld", p_ev[it*narr*nev+j*nev+iev]);
            }
#endif
            fprintf(fp, "\n");
        }
    }

    fclose(fp);
    free(p_ns);


#if defined(USE_LIKWID) || defined(USE_PAPIX6)
    free(p_ev);
#endif

#endif

    if (myrank == 0) {
        printf("Done. %f\n", y[narr/2][narr/2]);
    }

    for (size_t i = 0; i < narr; i ++) {
        free(x[i]);
    }
    for (size_t i = 0; i < narr; i ++) {
        free(y[i]);
    }

    free(x);
    free(y);


    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}
