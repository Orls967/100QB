#define _POSIX_C_SOURCE 200809L
#include "ising.h"
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static void build_temp_grid(double *T_out, int *n_out) {
    double raw[37];
    int k = 0;
    for (int i = 0; i < 25; ++i) {
        raw[k++] = 1.8 + (3.4 - 1.8) * (double)i / 24.0;
    }
    for (int i = 0; i < 12; ++i) {
        raw[k++] = 2.15 + (2.40 - 2.15) * (double)i / 11.0;
    }
    for (int i = 0; i < k - 1; ++i) {
        for (int j = i + 1; j < k; ++j) {
            if (raw[j] < raw[i]) {
                double t = raw[i];
                raw[i] = raw[j];
                raw[j] = t;
            }
        }
    }
    int n = 0;
    for (int i = 0; i < k; ++i) {
        if (n == 0 || fabs(raw[i] - T_out[n - 1]) > 1e-12) {
            T_out[n++] = raw[i];
        }
    }
    *n_out = n;
}

static void compute_stats_series(const double *x, int n, double *mean_out, double *err_out, double *tau_out) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) sum += x[i];
    double mean = sum / (double)n;
    double var = 0.0;
    for (int i = 0; i < n; ++i) {
        double d = x[i] - mean;
        var += d * d;
    }
    var /= (double)n;

    double tau = 0.5;
    if (var > 1e-30) {
        int max_lag = n / 4;
        if (max_lag > 2000) max_lag = 2000;
        for (int t = 1; t < max_lag; ++t) {
            double c = 0.0;
            for (int i = 0; i < n - t; ++i) {
                c += (x[i] - mean) * (x[i + t] - mean);
            }
            double rho = (c / (double)(n - t)) / var;
            if (rho <= 0.0) break;
            tau += rho;
            if ((double)t >= 6.0 * tau) break;
        }
    }
    double n_eff = (double)n / (2.0 * tau);
    if (n_eff < 1.0) n_eff = 1.0;
    *mean_out = mean;
    *err_out = sqrt(var / n_eff);
    *tau_out = tau;
}

static void run_parity(const char *out_dir) {
    const int L = 8;
    const double beta = 0.44;
    const int steps = 3000;
    const double p_table[5] = {1.0, 1.0, 1.0, exp(-4.0 * beta), exp(-8.0 * beta)};
    const double p_add = 1.0 - exp(-2.0 * beta);
    const char *names[3] = {"metropolis", "wolff", "sw"};

    char sum_path[512];
    snprintf(sum_path, sizeof(sum_path), "%s/parity_c_summary.csv", out_dir);
    FILE *fsum = fopen(sum_path, "w");
    fprintf(fsum, "algo,final_fnv1a_hex,final_E,final_M\n");

    for (int algo = 0; algo < 3; ++algo) {
        IsingLattice lat;
        ising_init(&lat, L);
        ising_set_cold(&lat);
        Xoshiro256pp rng;
        xoshiro_seed(&rng, compute_task_seed(L, 0, algo, 0));

        char path[512];
        snprintf(path, sizeof(path), "%s/parity_c_%s.csv", out_dir, names[algo]);
        FILE *fp = fopen(path, "w");
        fprintf(fp, "step,E_sum,M_sum\n");

        int64_t E = 0, M = 0;
        for (int s = 1; s <= steps; ++s) {
            if (algo == 0) ising_metropolis_sweep_table(&lat, &rng, p_table);
            else if (algo == 1) ising_wolff_step(&lat, &rng, p_add);
            else ising_sw_step(&lat, &rng, p_add);
            ising_energy_mag(&lat, &E, &M);
            fprintf(fp, "%d,%lld,%lld\n", s, (long long)E, (long long)M);
        }
        fclose(fp);
        uint64_t h = fnv1a_64(lat.spins, (size_t)lat.N);
        fprintf(fsum, "%s,%016llx,%lld,%lld\n", names[algo], (unsigned long long)h, (long long)E, (long long)M);
        ising_free(&lat);
    }
    fclose(fsum);
}

static void run_enum(const char *out_dir) {
    const double betas[5] = {0.2, 0.35, 0.4406867935097715, 0.5, 0.7};
    char path[512];
    snprintf(path, sizeof(path), "%s/enum_c.csv", out_dir);
    FILE *fp = fopen(path, "w");
    fprintf(fp, "L,beta,Z,e_mean,m2_mean\n");
    for (int L = 4; L <= 5; ++L) {
        for (int i = 0; i < 5; ++i) {
            EnumResult r = ising_enumerate_exact(L, betas[i]);
            fprintf(fp, "%d,%.17g,%.17g,%.17g,%.17g\n", L, r.beta, r.Z, r.e_mean, r.m2_mean);
        }
    }
    fclose(fp);
}

typedef struct {
    int L;
    int t_idx;
    double T;
    int algo;
    int n_therm;
    int n_meas;
    double e_mean, e_err, abs_m_mean, abs_m_err, m2_mean, m2_err, m4_mean, m4_err;
    double binder, chi, cv, tau_e, tau_m;
    uint64_t final_hash;
} SimTask;

typedef struct {
    SimTask *tasks;
    int n_tasks;
    int next_idx;
    pthread_mutex_t mtx;
} WorkQueue;

static void execute_sim_task(SimTask *tk) {
    IsingLattice lat;
    ising_init(&lat, tk->L);
    Xoshiro256pp rng;
    xoshiro_seed(&rng, compute_task_seed(tk->L, tk->t_idx, tk->algo, 0));
    ising_set_hot(&lat, &rng);

    const double beta = 1.0 / tk->T;
    const double p_table[5] = {1.0, 1.0, 1.0, exp(-4.0 * beta), exp(-8.0 * beta)};
    const double p_add = 1.0 - exp(-2.0 * beta);

    for (int s = 0; s < tk->n_therm; ++s) {
        if (tk->algo == 0) ising_metropolis_sweep_table(&lat, &rng, p_table);
        else if (tk->algo == 1) ising_wolff_step(&lat, &rng, p_add);
        else ising_sw_step(&lat, &rng, p_add);
    }

    const int M = tk->n_meas;
    double *es = (double *)malloc((size_t)M * sizeof(double));
    double *ams = (double *)malloc((size_t)M * sizeof(double));
    double *m2s = (double *)malloc((size_t)M * sizeof(double));
    double *m4s = (double *)malloc((size_t)M * sizeof(double));
    double sum_e2 = 0.0;
    const double invN = 1.0 / (double)lat.N;

    for (int s = 0; s < M; ++s) {
        if (tk->algo == 0) ising_metropolis_sweep_table(&lat, &rng, p_table);
        else if (tk->algo == 1) ising_wolff_step(&lat, &rng, p_add);
        else ising_sw_step(&lat, &rng, p_add);
        int64_t Ei = 0, Mi = 0;
        ising_energy_mag(&lat, &Ei, &Mi);
        double e = (double)Ei * invN;
        double m = (double)Mi * invN;
        double m2 = m * m;
        es[s] = e;
        ams[s] = fabs(m);
        m2s[s] = m2;
        m4s[s] = m2 * m2;
        sum_e2 += e * e;
    }

    compute_stats_series(es, M, &tk->e_mean, &tk->e_err, &tk->tau_e);
    compute_stats_series(ams, M, &tk->abs_m_mean, &tk->abs_m_err, &tk->tau_m);
    double dummy_t;
    compute_stats_series(m2s, M, &tk->m2_mean, &tk->m2_err, &dummy_t);
    compute_stats_series(m4s, M, &tk->m4_mean, &tk->m4_err, &dummy_t);

    tk->binder = 1.0 - tk->m4_mean / (3.0 * tk->m2_mean * tk->m2_mean);
    tk->chi = (double)lat.N * (tk->m2_mean - tk->abs_m_mean * tk->abs_m_mean) / tk->T;
    tk->cv = (double)lat.N * (sum_e2 / (double)M - tk->e_mean * tk->e_mean) / (tk->T * tk->T);
    tk->final_hash = fnv1a_64(lat.spins, (size_t)lat.N);

    free(es);
    free(ams);
    free(m2s);
    free(m4s);
    ising_free(&lat);
}

static void *worker_fn(void *arg) {
    WorkQueue *q = (WorkQueue *)arg;
    while (1) {
        pthread_mutex_lock(&q->mtx);
        int idx = q->next_idx++;
        pthread_mutex_unlock(&q->mtx);
        if (idx >= q->n_tasks) break;
        execute_sim_task(&q->tasks[idx]);
    }
    return NULL;
}

static void run_parallel_tasks(SimTask *tasks, int n_tasks, int n_threads) {
    if (n_threads <= 1) {
        for (int i = 0; i < n_tasks; ++i) execute_sim_task(&tasks[i]);
        return;
    }
    WorkQueue q;
    q.tasks = tasks;
    q.n_tasks = n_tasks;
    q.next_idx = 0;
    pthread_mutex_init(&q.mtx, NULL);
    pthread_t *tids = (pthread_t *)malloc((size_t)n_threads * sizeof(pthread_t));
    for (int t = 0; t < n_threads; ++t) pthread_create(&tids[t], NULL, worker_fn, &q);
    for (int t = 0; t < n_threads; ++t) pthread_join(tids[t], NULL);
    free(tids);
    pthread_mutex_destroy(&q.mtx);
}

static void run_threads_test(const char *out_dir, int n_threads) {
    double temps[4] = {2.0, 2.269185314213022, 2.5, 3.0};
    SimTask tasks[12];
    int k = 0;
    for (int ti = 0; ti < 4; ++ti) {
        for (int algo = 0; algo < 3; ++algo) {
            tasks[k].L = 16;
            tasks[k].t_idx = ti;
            tasks[k].T = temps[ti];
            tasks[k].algo = algo;
            tasks[k].n_therm = 500;
            tasks[k].n_meas = 2000;
            k++;
        }
    }
    run_parallel_tasks(tasks, k, n_threads);
    char path[512];
    snprintf(path, sizeof(path), "%s/threads_c_t%d.csv", out_dir, n_threads);
    FILE *fp = fopen(path, "w");
    fprintf(fp, "L,t_idx,algo,e_mean,m2_mean,final_hash\n");
    for (int i = 0; i < k; ++i) {
        fprintf(fp, "%d,%d,%d,%.17g,%.17g,%016llx\n", tasks[i].L, tasks[i].t_idx, tasks[i].algo,
                tasks[i].e_mean, tasks[i].m2_mean, (unsigned long long)tasks[i].final_hash);
    }
    fclose(fp);
}

static void run_full(const char *out_dir, int n_threads) {
    const int Ls[5] = {8, 16, 32, 64, 128};
    double temps[37];
    int n_temps = 0;
    build_temp_grid(temps, &n_temps);

    int max_tasks = 5 * n_temps * 3;
    SimTask *tasks = (SimTask *)calloc((size_t)max_tasks, sizeof(SimTask));
    int n_tasks = 0;
    for (int li = 0; li < 5; ++li) {
        int L = Ls[li];
        for (int ti = 0; ti < n_temps; ++ti) {
            for (int algo = 0; algo < 3; ++algo) {
                SimTask *tk = &tasks[n_tasks++];
                tk->L = L;
                tk->t_idx = ti;
                tk->T = temps[ti];
                tk->algo = algo;
                if (algo == 0) {
                    tk->n_therm = (L <= 64) ? 16384 : 8192;
                    tk->n_meas = (L <= 64) ? 262144 : 65536;
                } else if (algo == 1) {
                    tk->n_therm = 2048;
                    tk->n_meas = (L <= 64) ? 131072 : 32768;
                } else {
                    tk->n_therm = 2048;
                    tk->n_meas = (L <= 64) ? 131072 : 16384;
                }
            }
        }
    }

    double t0 = now_sec();
    run_parallel_tasks(tasks, n_tasks, n_threads);
    double dt = now_sec() - t0;
    printf("C full run (%d tasks, %d threads) finished in %.2f s\n", n_tasks, n_threads, dt);

    const char *names[3] = {"metropolis", "wolff", "sw"};
    char path[512];
    snprintf(path, sizeof(path), "%s/full_c.csv", out_dir);
    FILE *fp = fopen(path, "w");
    fprintf(fp, "L,t_idx,T,algo,n_therm,n_meas,e_mean,e_err,abs_m_mean,abs_m_err,m2_mean,m2_err,m4_mean,m4_err,binder,chi,cv,tau_e,tau_m,final_hash\n");
    for (int i = 0; i < n_tasks; ++i) {
        SimTask *tk = &tasks[i];
        fprintf(fp, "%d,%d,%.17g,%s,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%016llx\n",
                tk->L, tk->t_idx, tk->T, names[tk->algo], tk->n_therm, tk->n_meas,
                tk->e_mean, tk->e_err, tk->abs_m_mean, tk->abs_m_err,
                tk->m2_mean, tk->m2_err, tk->m4_mean, tk->m4_err,
                tk->binder, tk->chi, tk->cv, tk->tau_e, tk->tau_m, (unsigned long long)tk->final_hash);
    }
    fclose(fp);
    free(tasks);

    /* Save 3 L=128 snapshots at T=2.0, Tc, 3.0 */
    double snap_T[3] = {2.0, 2.269185314213022, 3.0};
    const char *snap_lbl[3] = {"T2.0", "Tc", "T3.0"};
    for (int si = 0; si < 3; ++si) {
        IsingLattice lat;
        ising_init(&lat, 128);
        Xoshiro256pp rng;
        xoshiro_seed(&rng, compute_task_seed(128, 900 + si, 1, 0));
        ising_set_hot(&lat, &rng);
        double p_add = 1.0 - exp(-2.0 / snap_T[si]);
        for (int s = 0; s < 1000; ++s) ising_wolff_step(&lat, &rng, p_add);
        snprintf(path, sizeof(path), "%s/snapshot_L128_%s.csv", out_dir, snap_lbl[si]);
        FILE *fs = fopen(path, "w");
        for (int y = 0; y < 128; ++y) {
            for (int x = 0; x < 128; ++x) {
                fprintf(fs, "%d%s", (int)lat.spins[y * 128 + x], (x + 1 == 128) ? "\n" : ",");
            }
        }
        fclose(fs);
        ising_free(&lat);
    }
}

static void run_bench(const char *out_dir) {
    const int L = 32;
    const double T = 2.269185314213022;
    const double beta = 1.0 / T;
    const double p_table[5] = {1.0, 1.0, 1.0, exp(-4.0 * beta), exp(-8.0 * beta)};
    const double p_add = 1.0 - exp(-2.0 * beta);
    const int reps = 7;

    char path[512];
    snprintf(path, sizeof(path), "%s/bench_c.csv", out_dir);
    FILE *fp = fopen(path, "w");
    fprintf(fp, "lang,variant,rep,ns_per_spin_update\n");

    for (int r = 0; r < reps; ++r) {
        IsingLattice lat;
        ising_init(&lat, L);
        Xoshiro256pp rng;
        xoshiro_seed(&rng, compute_task_seed(L, 0, 0, r));
        ising_set_hot(&lat, &rng);
        int sweeps = 2000;
        double t0 = now_sec();
        for (int s = 0; s < sweeps; ++s) ising_metropolis_sweep_table(&lat, &rng, p_table);
        double dt = now_sec() - t0;
        fprintf(fp, "c,metropolis_table,%d,%.17g\n", r, (dt * 1e9) / ((double)sweeps * lat.N));

        ising_set_hot(&lat, &rng);
        t0 = now_sec();
        for (int s = 0; s < sweeps; ++s) ising_metropolis_sweep_computed(&lat, &rng, p_table);
        dt = now_sec() - t0;
        fprintf(fp, "c,metropolis_computed,%d,%.17g\n", r, (dt * 1e9) / ((double)sweeps * lat.N));

        ising_set_hot(&lat, &rng);
        int w_steps = 2000;
        for (int s = 0; s < 200; ++s) ising_wolff_step(&lat, &rng, p_add);
        t0 = now_sec();
        for (int s = 0; s < w_steps; ++s) ising_wolff_step(&lat, &rng, p_add);
        dt = now_sec() - t0;
        fprintf(fp, "c,wolff,%d,%.17g\n", r, (dt * 1e9) / ((double)w_steps * lat.N));

        ising_set_hot(&lat, &rng);
        int sw_steps = 1000;
        t0 = now_sec();
        for (int s = 0; s < sw_steps; ++s) ising_sw_step(&lat, &rng, p_add);
        dt = now_sec() - t0;
        fprintf(fp, "c,sw,%d,%.17g\n", r, (dt * 1e9) / ((double)sw_steps * lat.N));
        ising_free(&lat);
    }

    /* Multi-spin coding sanity check vs scalar at L=16, T=2.5 */
    double ems = 0.0, m2ms = 0.0;
    double t0 = now_sec();
    ising_multispin_metropolis(16, 1.0 / 2.5, 4000, MASTER_SEED + 77, &ems, &m2ms);
    double dt_ms = now_sec() - t0;
    printf("C multispin L=16 T=2.5 (64 replicas x 4000 sweeps in %.3fs): e=%.6f, m2=%.6f\n", dt_ms, ems, m2ms);
    fclose(fp);
}

int main(int argc, char **argv) {
    const char *mode = "parity";
    const char *out_dir = "data";
    int n_threads = 8;
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--mode") == 0 && i + 1 < argc) mode = argv[++i];
        else if (strcmp(argv[i], "--out-dir") == 0 && i + 1 < argc) out_dir = argv[++i];
        else if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc) n_threads = atoi(argv[++i]);
    }
    if (strcmp(mode, "parity") == 0) run_parity(out_dir);
    else if (strcmp(mode, "enum") == 0) run_enum(out_dir);
    else if (strcmp(mode, "threads") == 0) run_threads_test(out_dir, n_threads);
    else if (strcmp(mode, "full") == 0) run_full(out_dir, n_threads);
    else if (strcmp(mode, "bench") == 0) run_bench(out_dir);
    return 0;
}
