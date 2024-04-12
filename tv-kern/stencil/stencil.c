#define _GNU_SOURCE
#define _ISOC11_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <stdlib.h>
#include <unistd.h>
#include <sched.h>
#include "mpi.h"
#ifdef PAPI
#include "papi.h"
#endif

// Warmup for 1000ms.
#ifndef NWARM
#define NWARM 2000
#endif

// Number of tests for each interval
#ifndef NTEST
#define NTEST 10
#endif

#ifndef NARR
#define NARR 100
#endif

#ifndef NPASS
#define NPASS 0 
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
main(void) {
    uint64_t ntest;
    /* Vars for Stencil */
    double w[NARR][NARR], Di[NARR][NARR], p[NARR][NARR], Kx[NARR][NARR], Ky[NARR][NARR];
    double rx, ry, pw;
    uint64_t narr;
    struct timespec tv;
    uint64_t volatile nsec_st, nsec_en; // For warmup
    int myrank, mycpu, nrank, errid;

    errid = MPI_Init(NULL, NULL);
    if (errid != MPI_SUCCESS) {
        printf("Faild to init MPI.\n");
        exit(1);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    mycpu = sched_getcpu();

    narr = NARR;
    ntest = NPASS + NTEST;



    for (uint64_t i = 0; i < narr; i ++) {
        fill_random(w[i], narr);
        fill_random(Di[i], narr);
        fill_random(p[i], narr);
        fill_random(Kx[i], narr);
        fill_random(Ky[i], narr);
    }
    fill_random(&rx, 1);
    fill_random(&ry, 1);
    

    if (myrank == 0) {
        printf("NTEST=%llu, NARR=%llu, rx=%f, ry=%f \n", ntest, narr, rx, ry);
    }

    // Warm up
    do {
        if (myrank == 0) {
            printf("Warming up for %d ms.\n", NWARM);
        }
        struct timespec tv;
        uint64_t volatile sec, nsec; // For warmup
        clock_gettime(CLOCK_MONOTONIC, &tv);
        sec = tv.tv_sec;
        nsec = tv.tv_nsec;
        nsec = sec * 1e9 + nsec + NWARM * 1e6;
        // tl_cg_calc_w stencil kernel
        while (tv.tv_sec * 1e9 + tv.tv_nsec < nsec) {
            for (uint64_t k = 0; k < narr; k ++) {
                for (uint64_t j = 0; j < narr; j ++) {
                    w[k][j] = Di[k][j] * p[k][j]  \
                              - ry * (Ky[k+1][j] * p[k+1][j] + Ky[k][j] * p[k-1][j]) \
                              - rx * (Kx[k][j+1] * p[k][j+1] + Kx[k][j] * p[k][j-1]);
                }
            }
            clock_gettime(CLOCK_MONOTONIC, &tv);
        }
    } while (0);

    if (myrank == 0) {
        printf("Start running.\n");
    }
    pw = 0;
#ifdef STIMING
    uint64_t *p_ns, ns0, ns1;
    char fname[] = "stiming_time.csv";
    char tot_time_fname[] = "stiming_tot_time.csv";
    p_ns =  (uint64_t *)malloc((ntest * narr) * sizeof(uint64_t));

#elif PAPI
    uint64_t *p_ns, ns0, ns1;
    char fname[] = "papi_time.csv";
    char tot_time_fname[] = "papi_tot_time.csv";
    int eventset = PAPI_NULL;
    p_ns =  (uint64_t *)malloc((ntest * narr) * sizeof(uint64_t));
    PAPI_library_init(PAPI_VER_CURRENT);
    PAPI_create_eventset(&eventset);
    PAPI_start(eventset);

#elif CGT
    uint64_t *p_ns, ns0, ns1;
    struct timespec tv1;
    char fname[] = "cgt_time.csv";
    char tot_time_fname[] = "cgt_tot_time.csv";
    p_ns =  (uint64_t *)malloc((ntest * narr) * sizeof(uint64_t));
#elif WTIME
    uint64_t *p_ns, ns0, ns1;
    char fname[] = "wtime_time.csv";
    char tot_time_fname[] = "wtime_tot_time.csv";
    p_ns =  (uint64_t *)malloc((ntest * narr) * sizeof(uint64_t));

#else
    char tot_time_fname[] = "clean_tot_time.csv";
#endif


    MPI_Barrier(MPI_COMM_WORLD);
    clock_gettime(CLOCK_MONOTONIC, &tv);
    nsec_st = tv.tv_sec * 1e9 + tv.tv_nsec;
    for (int i = 0; i < ntest; i ++) {
        // Timing.
        for (uint64_t k = 0; k < narr; k ++) {
#ifdef STIMING
            _read_ns (ns0);
            _mfence;
#elif PAPI
            ns0 = PAPI_get_real_nsec();
#elif CGT
            clock_gettime(CLOCK_MONOTONIC, &tv1);
            ns0 = tv1.tv_sec * 1e9 + tv1.tv_nsec;
#elif WTIME
            ns0 = (uint64_t)(MPI_Wtime() * 1e9);
#endif
            // tl_cg_calc_w stencil kernel
            for (uint64_t j = 0; j < narr; j ++) {
                w[k][j] = Di[k][j] * p[k][j]  \
                          - ry * (Ky[k+1][j] * p[k+1][j] + Ky[k][j] * p[k-1][j]) \
                          - rx * (Kx[k][j+1] * p[k][j+1] + Kx[k][j] * p[k][j-1]);
            }
#ifdef STIMING
            _read_ns (ns1);
            _mfence;
            p_ns[i*narr+k] = ns1 - ns0;
#elif PAPI
            ns1 = PAPI_get_real_nsec();
            p_ns[i*narr+k] = ns1 - ns0;
#elif CGT
            clock_gettime(CLOCK_MONOTONIC, &tv1);
            ns1 = tv1.tv_sec * 1e9 + tv1.tv_nsec;
            p_ns[i*narr+k] = ns1 - ns0;
#elif WTIME
            ns1 = (uint64_t)(MPI_Wtime() * 1e9);
            p_ns[i*narr+k] = ns1 - ns0;
#endif
            for (uint64_t j = 0; j < narr; j ++) {
                pw = pw + w[k][j] * p[k][j];
            }
        }
    MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (myrank == 0) {
        clock_gettime(CLOCK_MONOTONIC, &tv);
        nsec_en = tv.tv_sec * 1e9 + tv.tv_nsec;
        printf("Run time: %llu ns\n", nsec_en - nsec_st);
        FILE *fp_tot = fopen(tot_time_fname, "w");
        fprintf(fp_tot, "%llu\n", nsec_en - nsec_st);
        fclose(fp_tot);
    }

    // Gathering timing results.
    uint64_t *p_ns_all;
    p_ns_all = (uint64_t *)malloc(nrank*ntest*narr*sizeof(uint64_t));
    MPI_Gather(p_ns, ntest*narr, MPI_UINT64_T, p_ns_all, ntest*narr, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    if (myrank == 0) {
        FILE *fp = fopen(fname, "w");
        for (int irank = 0; irank < nrank; irank ++) {
            uint64_t ibias = irank * ntest * narr;
            for (uint64_t i = NPASS * narr; i < ntest * narr; i ++) {
                uint64_t iall = ibias + i;
                fprintf(fp, "%d,%llu\n", irank, p_ns_all[iall]);
            }
        }
        fclose(fp);
    }
    free(p_ns_all);
    MPI_Barrier(MPI_COMM_WORLD);


#ifdef PAPI
    PAPI_shutdown();
#endif


    MPI_Finalize();

    return 0;
}
