#define _GNU_SOURCE
#define _ISOC11_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <stdlib.h>
#include <unistd.h>
#include <sched.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "mpi.h"

#if defined(USE_PAPI) || defined(USE_PAPIX6)
#include "papi.h"

#elif defined(USE_LIKWID)
#include "likwid-marker.h"
#define NEV 20

#endif

// Warmup for 1000ms.
#ifndef NWARM
#define NWARM 1000
#endif

// Ignoring some tests at the beginning
#ifndef NPASS
#define NPASS 2
#endif

// Number of tests for each interval
#ifndef NTEST
#define NTEST 10
#endif

#ifndef TBASE
#define TBASE 100
#endif

// Uniform: V1: Number of DSub between two bins; V2: Number of bins.
// Normal: V1: sigma, standard error
// Pareto: V1: alpha, Pareto exponent 
#ifndef V1
#define V1 100
#endif

// Number of intervals (>=1)
#ifndef V2
#define V2 20
#endif

#ifndef FSIZE
#define FSIZE 0
#endif

#ifndef LCUT
#define LCUT 1-0x7fffffff
#endif

#ifndef HCUT
#define HCUT 0x7fffffff-1
#endif

// Prefetching the benchmarking instructions.
#ifndef NPRECALC
#define NPRECALC 1
#endif

#ifndef NSAMP
#define NSAMP 0
#endif

// Timing macros
#ifdef USE_TSC
#ifdef __x86_64__
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

#define _mfence asm volatile("lfence"   "\n\t":::);
#endif
#endif

#ifdef USE_TSC
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
#endif


/**
 * @brief A desginated one-line computing kernel.
 */
void
flush_cache(double *pf_a, double *pf_b, double *pf_c, uint64_t npf) {
#ifdef INIT
    for (uint64_t i = 0; i < npf; i ++) {
        pf_a[i] = i;
    }

#elif TRIAD
    for (uint64_t i = 0; i < npf; i ++) {
        pf_a[i] = 0.42 * pf_b[i] + pf_c[i];
    }

#elif SCALE
    for (uint64_t i = 0; i < npf; i ++) {
        pf_a[i] = 1.0001 * pf_b[i];
        pf_b[i] = 1.0001 * pf_a[i];
    }

#endif
    return;
}

/**
 * @brief Reading urandom file for random seeds.
 */
int 
get_urandom(uint64_t *x) {
    FILE *fp = fopen("/dev/urandom", "r");
    if (fp == NULL) {
        printf("Failed to open random file.\n");
        return 1;
    } 

    fread(x, sizeof(uint64_t), 1, fp); 

    fclose(fp);

    return 0;
}

/**
 * @brief Generate a shuffled list to loop through all ADD array's lengths.
 */
void
gen_walklist(uint64_t *len_list) {
#ifndef WALK_FILE
#ifdef UNIFORM
    for (int it = 0; it < NTEST; it ++) {
        for (int i = 0; i < V2; i ++) {
            len_list[NPASS+it*V2+i] = TBASE + i * V1;
        }
    }
    // Shuffle
    // Random seeds using nanosec timestamp
    struct timespec tv;
    uint64_t sec, nsec;
    clock_gettime(CLOCK_MONOTONIC, &tv);
    sec = tv.tv_sec;
    nsec = tv.tv_nsec;
    nsec = sec * 1e9 + nsec + NWARM * 1e6;
    srand(nsec); 
    for (int i = NPASS; i < NPASS + NTEST * V2; i ++){
        int r = (int) (((float)rand() / (float)RAND_MAX) * (float)(NTEST * V2));
        uint64_t temp = len_list[i];
        // Randomly swapping two vkern loop length
        len_list[i] = len_list[NPASS+r];
        len_list[NPASS+r] = temp;
    }
    for (int i = 0; i < NPASS; i ++) {
        len_list[i] = TBASE + V1 * V2;
    }
#else
    double v1 = V1;
    uint64_t seed;
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_ranlux389;
    r = gsl_rng_alloc(T);
    get_urandom(&seed);
    gsl_rng_set(r, seed);
    for (int i = 0; i < NPASS+NTEST; i++) {
        double x = 0;
#ifdef NORMAL
        do {
            x = 1.0 + gsl_ran_gaussian(r, v1);
        } while (x < LCUT || x > HCUT);
#elif PARETO
        do {
            x = gsl_ran_pareto(r, v1, 1.0);
        } while (x > HCUT);
#endif
        len_list[i] = (int)(x*TBASE);
    }
    gsl_rng_free(r);
#endif
#else
    // Reading walk list from file
    FILE *fp = fopen("walk.csv", "r");
    for (int it = 0; it < NTEST + NPASS; it ++) {
        fscanf(fp, "%llu\n", len_list+it);
    }
    fclose(fp);
#endif
}

int
main(int argc, char **argv) {
    int ntest;
    register uint64_t a, b;
    uint64_t *p_len, *p_ns;
    uint64_t register ns0, ns1;
    uint64_t tbase=TBASE, fsize = FSIZE, npf = 0; // npf: flush arr length
    double *pf_a, *pf_b, *pf_c;
    int myrank = 0, nrank = 1, errid = 0;

#ifdef UNIFORM
    uint64_t v1 = V1, v2 = V2;

#elif NORMAL
    double v1 = V1;

#elif PARETO
    double v1 = V1;
#endif

    errid = MPI_Init(NULL, NULL);
    if (errid != MPI_SUCCESS) {
        printf("Faild to init MPI.\n");
        exit(1);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


#ifdef UNIFORM
    ntest = NPASS + NTEST * V2;
#else
    ntest = NTEST + NPASS;
#endif

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
    p_ev = (int64_t *)malloc(ntest * nev * sizeof(int64_t));
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
    double ev_vals_0[NEV]={0}, ev_vals_1[NEV]={0}, time;
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
    p_ev = (int64_t *)malloc(ntest * nev * sizeof(int64_t));
    for (int iev = 0; iev < nev; iev ++) {
        ev_vals_0[iev] = 0;
        ev_vals_1[iev] = 0;
    }

#elif USE_TSC
    double tsc_freq = atof(argv[1]);

#endif

    p_len = (uint64_t *)malloc(ntest * sizeof(uint64_t));
    p_ns = (uint64_t *)malloc(ntest * sizeof(uint64_t));
    if (fsize) {
        npf = FSIZE / sizeof(double);
        pf_a = (double *)malloc(npf * sizeof(double));
        pf_b = (double *)malloc(npf * sizeof(double));
        pf_c = (double *)malloc(npf * sizeof(double));
    }
    for (uint64_t i = 0; i < npf; i ++) {
        pf_a[i] = 1.1;
        pf_b[i] = 1.1;
        pf_c[i] = 1.1;
    }

    if (myrank == 0) {
#ifdef UNIFORM
        printf("Generating uniform ditribution. Tbase=%lu, interval=%lu, nint=%lu, Ntest=%u.\n",
            tbase, v1, v2, NTEST);

#elif NORMAL
        printf("Generating normal ditribution. Tbase=%lu, sigma=%.4f, Ntest=%lu.\n",
            tbase, v1, NTEST);

#elif PARETO
        printf("Generating Pareto ditribution. Tbase=%lu, alpha=%.4f, Ntest=%lu.\n",
            tbase, v1, NTEST);
#else
        printf("Reading vkern walk list from walk.csv.\n");

#endif
        fflush(stdout);
        gen_walklist(p_len);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(p_len, ntest, MPI_UINT64_T, 0, MPI_COMM_WORLD);


    a = 0;
    if (argc < 10) {
        b = 1;
    } else {
        b = atoi(argv[1]);
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

        while (tv.tv_sec * 1e9 + tv.tv_nsec < nsec) {
            a = b + 1;
            clock_gettime(CLOCK_MONOTONIC, &tv);
        }
    } while (0);

    if (myrank == 0) {
        printf("Start random walking.\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);
    ns1 = 0;
    ns0 = 0;
    MPI_Barrier(MPI_COMM_WORLD);
    for (int iwalk = 0; iwalk < ntest; iwalk ++) {
        register uint64_t n = p_len[iwalk];
        register uint64_t npre = NPRECALC;
        register uint64_t ra, rb, rc;
        // struct timespec tv;

        // Flushing
        if (npf) {
            flush_cache(pf_a, pf_b, pf_c, npf);
        }

        // Instruction preload.
        ra = n + npre;
        rb = b;
        while (ra != n) {
            ra = ra - rb;
        }

        // Timing.
        while (ra!= 0) {
            ra = ra - rb;
        }

        ra = NSAMP;
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
        ns0 = ns1;
        LIKWID_MARKER_START("vkern"); 
#elif USE_TSC
        tsc_start(&ns0);
#else
        _read_ns (ns0);
        _mfence;

#endif

        while (ra!= 0) {
            ra = ra - rb;
        }
#ifdef USE_PAPI
        ns1 = PAPI_get_real_nsec();

#elif USE_PAPIX6
        ns1 = PAPI_get_real_nsec();
        PAPI_read(eventset, ev_vals_1);
        for (int iev = 0; iev < nev; iev ++) {
            p_ev[iwalk * nev + iev] = (int64_t)(ev_vals_1[iev] - ev_vals_0[iev]);
        }

#elif USE_CGT
        clock_gettime(CLOCK_MONOTONIC, &tv);
        ns1 = tv.tv_sec * 1e9 + tv.tv_nsec;

#elif USE_WTIME
        ns1 = (uint64_t)(MPI_Wtime() * 1e9);

#elif USE_LIKWID
        LIKWID_MARKER_STOP("vkern"); 
        LIKWID_MARKER_GET("vkern", &nev, (double*)ev_vals_1, &time, &count);
        for (int iev = 0; iev < nev; iev ++) {
            p_ev[iwalk * nev + iev] = (int64_t)ev_vals_1[iev] - (int64_t)ev_vals_0[iev];
            ev_vals_0[iev] = ev_vals_1[iev];
        }
        ns1 = (uint64_t)(1e9*time);
#elif USE_TSC
        tsc_stop(&ns1);
        ns1 = (uint64_t)((double)(ns1 - ns0) / tsc_freq);
        ns0 = 0;
#else
        _read_ns (ns1);
        _mfence;
#endif
        p_ns[iwalk] = ns1 - ns0;
        
        a += ra;
        MPI_Barrier(MPI_COMM_WORLD);
    }

#ifdef USE_PAPI
    PAPI_shutdown();

#elif USE_PAPIX6
    PAPI_shutdown();

#elif USE_LIKWID
    // Finalize LIKWID
    LIKWID_MARKER_CLOSE;    
#endif
    MPI_Barrier(MPI_COMM_WORLD);

    // Each rank writes its own file.
    char fname[4096], myhost[1024];
    gethostname(myhost, 1024);
#ifdef USE_PAPI
    sprintf(fname, "papi_time_%d_%s_tf.csv", myrank, myhost);
    FILE *fp = fopen(fname, "w");
#elif USE_CGT
    sprintf(fname, "cgt_time_%d_%s_tf.csv", myrank, myhost);
    FILE *fp = fopen(fname, "w");
#elif USE_WTIME
    sprintf(fname, "wtime_time_%d_%s_tf.csv", myrank, myhost);
    FILE *fp = fopen(fname, "w");
#elif USE_PAPIX6
    sprintf(fname, "papix6_time_%d_%s_tf.csv", myrank, myhost);
    FILE *fp = fopen(fname, "w");
#elif USE_LIKWID
    sprintf(fname, "likwid_time_%d_%s_tf.csv", myrank, myhost);
    FILE *fp = fopen(fname, "w");
#elif USE_TSC
    sprintf(fname, "tsc_time_%d_%s.csv", myrank, myhost);
    FILE *fp = fopen(fname, "w");
#else
    sprintf(fname, "stiming_time_%d_%s_tf.csv", myrank, myhost);
    FILE *fp = fopen(fname, "w");
#endif
    for (int iwalk = NPASS; iwalk < ntest; iwalk ++) {
        fprintf(fp, "%d,%lu,%lu", myrank, NSAMP, p_ns[iwalk]);
#if defined(USE_LIKWID) || defined(USE_PAPIX6)
        for (int iev = 0; iev < nev; iev ++) {
            fprintf(fp, ",%ld", p_ev[iwalk*nev+iev]);
        }
#endif
        fprintf(fp, "\n");
    }
    fclose(fp);
    MPI_Barrier(MPI_COMM_WORLD);
#if defined(USE_LIKWID) || defined(USE_PAPIX6)
    free(p_ev);
#endif

    if (myrank == 0) {
        printf("Done. %lu\n", a);
    }

    MPI_Finalize();

    if (npf) {
        free(pf_a);
        free(pf_b);
        free(pf_c);
    }

    free(p_len);
    free(p_ns);

    return 0;
}
