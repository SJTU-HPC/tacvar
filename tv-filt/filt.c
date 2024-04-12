/**
 * @file filt.c
 * @author Key Liao
 *
 */
#define _GNU_SOURCE
#define _ISOC11_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#ifdef __linux__
#include <error.h>
#endif
#ifdef __APPLE__
#include <errno.h>
#endif
#include <math.h>
#include <sys/time.h>
#include "argp.h"
#include "filt.h"

#ifndef NTILE
#define NTILE 1000
#endif


/*=== ARGP ===*/
static struct argp_option options[] = {
    {"met-file", 'm', "FILE", 0, "Input measurement file", 0},
    {"sample-file", 's', "FILE", 0, "Input timing fluctuation file", 0},
    {"width", 'w', "TIME", 0, "The least interval of a time bin (ns).", 0},
    {"nsamp", 'n', "NSAMP", 0, "Number of samples in each optimization step", 0},
    {"plow", 'l', "PROB", 0, "Lowest threshold of possibility of a data bin", 0},
    {"cut-x", 'x', "PROB", 0, "Cut the highest probability of met array", 0},
    {"cut-y", 'y', "PROB", 0, "Cut the highest probability of timing fluctuation array", 0},
    {"cut-z", 'z', "PROB", 0, "Cut the highest probability in w-distance calculation", 0},
    {0}
};

// Parser function
static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    filt_param_t *args = state->input;

    switch (key) {
        case 'm':
            args->in_tm_file = arg;
            break;
        case 's':
            args->in_tf_file = arg;
            break;
        case 'w':
            args->width = atoi(arg);
            break;
        case 'n':
            args->nsamp = atoi(arg);
            break;
        case 'l':
            args->p_low = atof(arg);
            break;
        case 'x':
            args->p_xcut = atof(arg);
            break;
        case 'y':
            args->p_ycut = atof(arg);
            break;
        case 'z':
            args->p_zcut = atof(arg);
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

// argp parser structure
static struct argp argp = {options, parse_opt, NULL, NULL};


/*=== BEGIN: Implementations ===*/

int
cmp(const void *a, const void *b) {
    return (*(int*)a - *(int*)b);
}


int
sim_met(prob_hist_t *tmh, prob_hist_t *trh, i64 *tf_arr, u64 tf_len, double *psim, u64 nsamp) {
    u64 nbin = tmh->nbin;
    double *ptrhc, trpmax; // Cumulative probabilities of trh.
    struct timeval tv;

    ptrhc = (double *)malloc(nbin*sizeof(double));

    gettimeofday(&tv, NULL);
    srand(tv.tv_sec * 1e6 + tv.tv_usec);

    // Calculating tr's cumulative probability
    ptrhc[0] = 0.0;
    psim[0] = 0;
    for (u64 ir = 1; ir < nbin; ir ++) {
        ptrhc[ir] = ptrhc[ir-1] + trh->pbin[ir-1].p;
        psim[ir] = 0;
    }
    trpmax = ptrhc[nbin-1];
    
    for (u64 isamp = 0; isamp < nsamp; isamp ++){
        register double px = (double)rand() / (double)RAND_MAX;
        if (px < trpmax) {
            register i64 tf = tf_arr[(u64)((double)rand() / (double)RAND_MAX * (double)tf_len)];
            for (u64 ir = 0; ir < nbin - 1; ir ++) {
                if (px < ptrhc[ir+1]) {
                    register i64 tsim;
                    register double p0, p1, t0, t1, tx;
                    p0 = ptrhc[ir];
                    p1 = ptrhc[ir+1];
                    t0 = trh->pbin[ir].t;
                    t1 = trh->pbin[ir+1].t;
                    tx = t0 + (px - p0) * (t1 - t0) / (p1 - p0);
                    tsim = (i64)tx + tf;
                    if (tsim >= trh->pbin[nbin-1].t) {
                        psim[nbin-1] ++;
                    } else {
                        for (u64 im = 0; im < nbin - 1; im ++) {
                            if (tsim < tmh->pbin[im+1].t) {
                                psim[im] ++;
                                break;
                            }
                        }
                    }
                    break;
                }
            }
        }
    }

    for (u64 isim = 0; isim < nbin; isim ++) {
        // printf("%f ", psim[isim]);
        psim[isim] = psim[isim] / (double)nsamp;
    }

    free(ptrhc);

    return 0;
}


int
sim_verify(prob_hist_t *tmh, prob_hist_t *trh, i64 *tf_arr, u64 tf_len, i64 *sim_cdf, u64 nsamp) {
    u64 nbin = tmh->nbin;
    i64 *sim_arr;
    double *ptrhc, trpmax; // Cumulative probabilities of trh.
    struct timeval tv;

    ptrhc = (double *)malloc(nbin*sizeof(double));
    sim_arr = (i64 *)malloc(nsamp*sizeof(i64));

    gettimeofday(&tv, NULL);
    srand(tv.tv_sec * 1e6 + tv.tv_usec);

    // Calculating tr's cumulative probability
    ptrhc[0] = 0.0;
    for (u64 ir = 1; ir < nbin; ir ++) {
        ptrhc[ir] = ptrhc[ir-1] + trh->pbin[ir-1].p;
    }
    trpmax = ptrhc[nbin-1];
    
    for (u64 isamp = 0; isamp < nsamp; isamp ++){
        // mul trpmax instead of 1, to tackle the little error from trpmax to 1.0
        register double px = (double)rand() / (double)RAND_MAX * trpmax;
        register i64 tf = tf_arr[(u64)((double)rand() / (double)RAND_MAX * (double)tf_len)];
        for (u64 ir = 0; ir < nbin - 1; ir ++) {
            if (px < ptrhc[ir+1]) {
                register i64 tsim;
                register double p0, p1, t0, t1, tx;
                p0 = ptrhc[ir];
                p1 = ptrhc[ir+1];
                t0 = trh->pbin[ir].t;
                t1 = trh->pbin[ir+1].t;
                tx = t0 + (px - p0) * (t1 - t0) / (p1 - p0);
                sim_arr[isamp] = (i64)tx + tf;
                break;
            }
        }
    }

    qsort(sim_arr, nsamp, sizeof(i64), cmp);

    for (int i = 0; i < NTILE; i ++) {
        sim_cdf[i] = sim_arr[abs((int)(((double)i/(double)NTILE)*(double)nsamp))];
    }
    sim_cdf[NTILE] = sim_arr[nsamp-1];

    free(sim_arr);
    free(ptrhc);

    return 0;
}

int
calc_tr(prob_hist_t *tmh, prob_hist_t *trh, i64 *tf_arr, u64 tf_len, filt_param_t *args){
		
	u64 nbin, imb;
	int err = 0;
	double tot_p = 0.0;	// Tracking total probability in trh for the exit condition
	double simh[tmh->nbin], simh2[tmh->nbin];	// Simulated probabilities
    double dp_min;

    dp_min = fabs(0.01 * args->p_low); // The min probability gap of optimization

	// Init trh
	nbin = tmh->nbin;
	trh->nbin = nbin;
	trh->pbin = (prob_bin_t *)malloc(nbin * sizeof(prob_bin_t));
	if (trh->pbin == NULL) {
        printf("[FilT-calc_tr] ERROR, trh->pbin allocation failed.\n");
        err = errno;
        return err;
    }
	for (u64 i = 0; i < nbin; i ++) {
		trh->pbin[i].t = tmh->pbin[i].t - tf_arr[0];
		trh->pbin[i].p = 0.0;
        simh[i] = 0;
	}

	// Loop over each time section in trh
    imb = 0;
	while (imb < nbin - 1) {
		i64 tr_l = trh->pbin[imb].t, tr_r = trh->pbin[imb+1].t; // tr in [tr_l, tr_r)
        double p, tr_p = 0;
        double dp, dmdr, dp0, dp1;

        // ====== S1: Initial condition ======
        p = tmh->pbin[imb].p - simh[imb];
        if (p < dp_min) {
            trh->pbin[imb].p = 0;
            imb += 1;
            continue;
        }
        trh->pbin[imb].p = p;
        printf("[FilT-calc_tr] IMB=%lld, [%lld, %lld), tr_p=%f\n", imb, tr_l, tr_r, p);
        fflush(stdout);

        // ====== S2: Optimization ======
        sim_met(tmh, trh, tf_arr, tf_len, simh2, args->nsamp);
        dp = simh2[imb] - tmh->pbin[imb].p;
        dmdr = (simh2[imb] - simh[imb]) / trh->pbin[imb].p;
        printf("[FilT-calc_tr] bin_p=%f, delta_p=%f, dmdr=%f\n", trh->pbin[imb].p, dp, dmdr);
        dp0 = 0x7fffffff;
        dp1 = dp;
        do {
            tr_p = trh->pbin[imb].p;
            trh->pbin[imb].p -= dp1 / dmdr;
            dp0 = dp1;
            for (u64 isimh = 0; isimh < nbin; isimh ++) {
                simh[isimh] = simh2[isimh];
            }
            if (trh->pbin[imb].p <= 0) {
                // We do not want nagative probability
                dp1 = 0; // Don't sim, exit.
                break;
            } else {
                sim_met(tmh, trh, tf_arr, tf_len, simh2, args->nsamp);
                dp1 = simh2[imb] - tmh->pbin[imb].p;
                printf("[FilT-calc_tr] bin_p=%f, delta_p=%f\n", trh->pbin[imb].p, dp1);
            }
            // If the gap <= dp_min, no further optimization.
            if (fabs(dp1) <= dp_min) {
                tr_p = trh->pbin[imb].p;
                break;
            }
        } while (fabs(dp0) > fabs(dp1));
        tot_p += tr_p;
        // exit condition
        if (tot_p >= 1.0) {
            // If the last try is closer to 1, then using the last tr_p as the final bin.
            if (fabs(tot_p - 1) > fabs(tot_p - tr_p - 1)) {
                tot_p -= tr_p;
                trh->pbin[imb].p = 0;
            } else {
                trh->pbin[imb].p = tr_p;
            }
            break;
        } else {
            trh->pbin[imb].p = tr_p;
            printf("[FilT-calc_tr] p=%f, tot_p=%f\n", trh->pbin[imb].p, tot_p);
            fflush(stdout);
            imb ++;
        }
	}
    ep = fabs(tot_p - 1);
    printf("[FilT-calc_tr] tot_p=%f, Normalizing probabilities...", tot_p);
    fflush(stdout);
    for (u64 i = 0; i < trh->nbin; i ++) {
        trh->pbin[i].p /= tot_p;
    } 
    printf("Done.\n");
    fflush(stdout);

	return err;
}

int
slice(i64 *arr, u64 len, double p_low, i64 width, prob_hist_t *phist){
    struct timeval st, en;
    int err = 0;
    i64 tmin, tmax; // The max and min in arr.
    u64 ibin, nbin, nc_low, nbv, nc; // init bin count, lowest point count, valid nbin
    i64 *pbins;     // Left edge of each bin.
    u64 *pcnts;    // Point counts of each bin.

    // For a phist = [[t0, p0], [t1, p1], ..., [t(n-1), p(n-1)]], it means for a
    // random measurement result t, the probability t in [t0, t1) is p0.
    // For any p in phist, p >= p_low except for p(n-1)≡0
    // t[i] - t[i-1] is always a positive multiple to width

    // Counting
    tmin = arr[0];
    tmax = arr[len-1];
    
    nbin = (u64)((float)(tmax - tmin) / (float)width);
    nbin = nbin < 2? 2: nbin + 1;
    printf("[FilT-slice] nbin=%llu, tmin=%lld, tmax=%lld \n", nbin, tmin, tmax);

    // Corner case, only 1 or 2 bin
    if (nbin <= 2) {
        phist->nbin = nbin;
        phist->pbin = (prob_bin_t *)malloc(nbin * sizeof(prob_bin_t));
        if (phist->pbin == NULL) {
            printf("[FilT-slice] ERROR, phist->pbin allocation failed.\n");
            err = errno;
            return err;
        }
		phist->pbin[nbin-1].t = tmax + 1;
        phist->pbin[nbin-1].p = 0;
        phist->pbin[0].t = tmin;
        phist->pbin[0].p = 1.0;

		return err;
    }
    
    pbins = (i64 *)malloc(nbin * sizeof(i64));
    pcnts = (u64 *)malloc(nbin * sizeof(u64));
    if (pbins == NULL || pcnts == NULL) {
        printf("[FilT-slice] ERROR, pbins and pcnts allocation failed.\n");
        err = errno;
        return err;
    }
    // Initalizing hist
    printf("[FilT-slice] Binning ... ");
    for (u64 i = 0; i < nbin; i ++) {
        pbins[i] = tmin + i * width;
        pcnts[i] = 0;
    }

    // Counting in even bins
    ibin = 0;   // current bin id
    pcnts[0] = 1;
    for (u64 i = 1; i < len; i ++) {
        if (arr[i] == arr[i-1]) {
            pcnts[ibin] ++;
        } else {
            ibin = (u64)((float)(arr[i] - tmin) / (float)width);
            pcnts[ibin] ++;
        }
    }

    // Merging for p_low
    nc_low = (u64)((double)len * p_low) + 1;
    nbv = 0;
    nc = 0;
	// The last bin is always empty, only for marking the right edge of the histo.
    for (u64 i = 0; i < nbin-1; i ++) {
        nc += pcnts[i];
        if (nc >= nc_low) {
            pcnts[nbv] = nc;
            nc = 0;
            nbv ++;
            pbins[nbv] = pbins[i+1];
        }
    }
    pcnts[nbv-1] += nc;
    pbins[nbv] = pbins[nbin-1];
    pcnts[nbv] = 0;
    nbv ++;

    // Calculating probability
    phist->nbin = nbv;
    phist->pbin = (prob_bin_t *)malloc(nbv * sizeof(prob_bin_t));
    if (phist->pbin == NULL) {
        printf("[FilT-slice] ERROR, phist->pbin allocation failed.\n");
        err = errno;
        return err;
    }
    for (u64 i = 0; i < nbv; i ++) {
        phist->pbin[i].t = pbins[i];
        phist->pbin[i].p = (double)pcnts[i] / (double)len;
    }
    printf("Done.\n");
    printf("[FilT-slice] Final bin count: %llu\n", phist->nbin);

    free(pbins);
    free(pcnts);

    return err;
}

void
calc_w(i64 *tm_arr, u64 tm_len, i64 *sim_cdf, i64 *w_arr, double *wp_arr, double p_zcut) {
    double w = 0, wm = 0;
    int wtile = (int)((1 - p_zcut) * (double)NTILE);
    for (size_t i = 0; i < NTILE; i ++) {
        size_t itm = (size_t)((double)i / (double)NTILE * tm_len);
        w_arr[i] = sim_cdf[i] - tm_arr[itm];
        wp_arr[i] = (double)w_arr[i] / tm_arr[itm];
    }
    for (size_t i = 0; i < wtile; i ++) {
        size_t itm = (size_t)((double)i / (double)NTILE * tm_len);
        w += abs(w_arr[i]);
        wm += tm_arr[itm];
    }
    w = w / (double)wtile;
    wm = wm / (double)wtile;
    er = w / wm;

    printf(" W-Distance=%f  ", w);
    FILE *fp=fopen("wd.out", "w");
    fprintf(fp, "%f", w);

    return;
}

int
read_csv(char *fpath, double pcut, u64 *len, i64 **arr) {
    FILE *fp;
    i64 fsize;
    i64 nrow;
    char *buf;
    int err = 0;
    i64 i = 0;
    i64 irow = 0, ist = 0;

    printf("[FilT-parse_csv] Reading %s.\n", fpath);
    fp = fopen(fpath, "r");
    if (fp == NULL) {
        printf("[FilT-parse_csv] ERROR, file does not exist at given path: %s.\n", fpath);
        return errno;
    }
    fseek(fp, 0, SEEK_END);
    fsize = ftell(fp);
    rewind(fp);

    // Allocate memory for raw file.
    buf = (char *)malloc(fsize + 1);
    if (buf == NULL) {
        printf("[FilT-parse_csv] ERROR, read buffer allocation failed.\n");
        fclose(fp);
        return errno;
    }

    // Read the file content into the character array
    fread(buf, fsize, 1, fp);
    buf[fsize] = '\0';  // null terminator
    fclose(fp);

    // Counting how many rows in the file.
    i = 0;
    irow = 0;
    while (buf[i]!='\0') {
        if (buf[i] == '\n'){
            irow ++;
        }
        i ++;
    }
    nrow = irow;
    *len = (u64)((double)nrow * (1.0 - pcut));
    printf("[FilT-parse_csv] %lld data points, drop the largest %f.\n", nrow, pcut);

    // Allocating memory space for raw data array.
    *arr = (i64 *)malloc(nrow*sizeof(i64));
    if (*arr == NULL) {
        printf("[FilT-parse_csv] ERROR, read buffer allocation failed.\n");
        fclose(fp);
        err = errno;
        return err;
    }
    // Converting char to i64;
    i = 0;
    irow = 0;
    ist = 0;
    while (buf[i]!='\0') {
        // Loop over lines.
        if (buf[i] == '\n') {
            buf[i] = '\0';
            (*arr)[irow] = strtoll(&(buf[ist]), NULL, 10);
            irow ++;
            ist = i + 1;
        }
        i ++;
    }
    qsort(*arr, nrow, sizeof(i64), cmp);

    return 0;
}

int
run_filt(filt_param_t *args, prob_hist_t *tr_hist, i64 *sim_cdf, 
            i64 *w_arr, double *wp_arr){
    u64 tm_len, tf_len;
    i64 *tm_arr, *tf_arr;
    prob_hist_t tm_hist;
    int err;

    // Parsing csv files and slicing specified column into histogram, saving to pmet_hist and ptf_hist
    // malloc inside slice function
    printf("[FilT-run_filt] Parsing measurement file %s\n", args->in_tm_file);
    err = read_csv(args->in_tm_file, args->p_xcut, &tm_len, &tm_arr);
    if (err) {
        printf("[FilT-run_filt] Error in reading measurement csv file. ERRCODE %d\n", err);
        return err;
    }

    printf("[FilT-run_filt] Slicing measurement file %s\n", args->in_tm_file);
    err = slice(tm_arr, tm_len, args->p_low, args->width, &tm_hist);
    if (err) {
        printf("[FilT-run_filt] Error in slicing met array. ERRCODE %d\n", err);
        return err;
    }
    // Print measured hist for debugging.
    printf("[FilT-run_filt] Histogram of measured run times:\ntime\t\tp\n");
    for (size_t i = 0; i < tm_hist.nbin - 1; i ++) {
        printf("[%ld, %ld)\t%.7f\n", 
                tm_hist.pbin[i].t, tm_hist.pbin[i+1].t, tm_hist.pbin[i].p);
    }

    printf("[FilT-run_filt] Parsing timing fluctuation file %s\n", args->in_tf_file);
    err = read_csv(args->in_tf_file, args->p_ycut, &tf_len, &tf_arr);
    if (err) {
        printf("[FilT-run_filt] Error in parsing timing fluctuations csv file. ERRCODE %d\n", err);
        return err;
    }

    printf("[FilT-run_filt] Start estimating real run time distribution.\n");
    err = calc_tr(&tm_hist, tr_hist, tf_arr, tf_len, args);
    if (err) {
        printf("[FilT-run_filt] Error in transposed convolution. ERRCODE %d\n", err);
        return err;
    }

    printf("[FilT-run_filt] Verifying the estimation...");
    sim_verify(&tm_hist, tr_hist, tf_arr, tf_len, sim_cdf, args->nsamp);
    calc_w(tm_arr, tm_len, sim_cdf, w_arr, wp_arr, args->p_zcut);
    printf("Done.\n");

    free(tm_arr);
    free(tf_arr);
    free(tm_hist.pbin);
    return 0;
}

//int
//arg_parse(int argc, char **argv){
//
//}

/*=== END: Implementations ===*/

/*=== BEGIN: Main Entry ===*/

int
main(int argc, char **argv){
    filt_param_t args;
    prob_hist_t tr_hist;
    int err;
    i64 sim_cdf[NTILE+1], w_arr[NTILE+1];
    double wp_arr[NTILE+1];
    
    args.in_tm_file = "met.csv";
    args.in_tf_file = "tf.csv";
    args.out_trh_file = "tr_hist.csv";
    args.out_sim_file = "sim_cdf.csv";
    args.width = 100;
    args.nsamp = 1000000;
    args.p_low = 0.001;
    args.p_xcut = 0.0;
    args.p_ycut = 0.0;

    // Parse command line
    argp_parse(&argp, argc, argv, 0, 0, &args);

    printf("[FilT-main] Met in file: %s\n", args.in_tm_file);
    printf("[FilT-main] Timing fluctuation file: %s\n", args.in_tf_file);
    printf("[FilT-main] p_low=%f, width=%llu, nsamp=%llu, "
            "met_cut=%f, tf_cut=%f\n", 
            args.p_low, args.width, args.nsamp, args.p_xcut, args.p_ycut);

    // Reading input files and estimating the real run time distribution.
    err = run_filt(&args, &tr_hist, sim_cdf, w_arr, wp_arr);

    if (err == 0) {
        printf("[FilT-main] Writing results to %s...", args.out_trh_file);
        FILE *fp = fopen(args.out_trh_file, "w");
        for (i64 i = 0; i < tr_hist.nbin; i ++) {
            fprintf(fp, "%lld, %lf\n", tr_hist.pbin[i].t, tr_hist.pbin[i].p);
        }
        fclose(fp);
        free(tr_hist.pbin);
        printf("Done.\n");
        printf("[FilT-main] Writing verification simulation results to %s...", args.out_sim_file);
        fp = fopen(args.out_sim_file, "w");
        for (int i = 0; i < NTILE+1; i ++) {
            fprintf(fp, "%d, %ld, %ld, %lf\n", i, sim_cdf[i], w_arr[i], wp_arr[i]);
        }
        fclose(fp);
        fp = fopen("er.out", "w");
        fprintf(fp, "%f", er);
        fclose(fp);
        fp = fopen("ep.out", "w");
        fprintf(fp, "%f", ep);
        fclose(fp);
        printf("Done. er=%f, ep=%f\n", er, ep);
    } else {
        printf("[FilT-main] run_filt returned with errors. ERRCODE %d\n", err);
        return err;
    }

    return 0;
}

/*=== END: Main Entry ===*/
