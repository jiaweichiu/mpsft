#include <fstream>
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <glog/logging.h>

#include "computefourier.h"
#include "fft.h"
#include "filters.h"
#include "parameters.h"
#include "plot.h"
#include "timer.h"
#include "utils.h"

// JW: our globals
// gGraphType=1: Fix S, increase N
// gGraphType=2: Fix N, increase S
// gGraphType=3: Fix S,N, increase SNR.
// gVer=1: sFFT1.0;  gVer=2: sFFT2.0
int gOuterReps, gInnerReps, gComputeError, gGraphType, gVer;
/*  gOuterReps=atoi(argv[1]);
  gInnerReps=atoi(argv[2]);
  gComputeError=atoi(argv[3]);
  gGraphType=atoi(argv[4]);
  gVer=atoi(argv[5]);*/

/*
  Run an experiment.  Parameters are:

  x: the signal
  n: the length of x

  lobefrac_loc:   during location, the main lobe of the filter has half
                  width n*lobefrac_loc
  tolerance_loc:  the linf norm of the residuals of the filter
  b_loc:          the number of adjacent filters to add

  B_loc:          number of samples in subsampling during location
  B_thresh:       number of samples considered "heavy"
  loops_loc:      number of location loops
  loops_thresh:   number of times a coordinate must seem heavy.

  *_est:          ditto as above, but for estimation.

  repetitions:    repeat the experiment this many times for timing
  LARGE_FREQ:     locations of largest coefficients.
  k:              number of HHs, used for evaluation only.
  x_f:            true fft.
 */
void run_experiment(complex_t *x, int n, double lobefrac_loc,
                    double tolerance_loc, int b_loc, int B_loc, int B_thresh,
                    int loops_loc, int loops_thresh, double lobefrac_est,
                    double tolerance_est, int b_est, int B_est, int loops_est,
                    int W_Comb, int Comb_loops, int repetitions, bool FFTW_OPT,
                    int *LARGE_FREQ, int k, complex_t *x_f, double &SFFT_time,
                    double &FFTW_time, double &error_per_entry) {
  if (WITH_COMB)
    assert(B_thresh < W_Comb);
  if (ALGORITHM1) {
    assert(B_thresh < B_loc);
    assert(loops_thresh <= loops_loc);
  }

  int w_loc;
  complex_t *filtert =
      make_dolphchebyshev_t(lobefrac_loc, tolerance_loc, w_loc);
  Filter filter = make_multiple_t(filtert, w_loc, n, b_loc);

  int w_est;
  complex_t *filtert_est =
      make_dolphchebyshev_t(lobefrac_est, tolerance_est, w_est);
  Filter filter_est = make_multiple_t(filtert_est, w_est, n, b_est);

  real_t filter_noise = 0, filter_noise_est = 0;
  for (int i = 0; i < 10; i++) {
    filter_noise =
        std::max(filter_noise, std::max(cabs(filter.freq[n / 2 + i]),
                                        cabs(filter.freq[n / 2 - i])));
    filter_noise_est =
        std::max(filter_noise_est, std::max(cabs(filter_est.freq[n / 2 + i]),
                                            cabs(filter_est.freq[n / 2 - i])));
  }

#ifdef DEBUG
  complex_t *tmp = (complex_t *)calloc(n, sizeof(*tmp));
  memcpy(tmp, filter.time, w_loc * sizeof(*tmp));
  free(filter.time);
  filter.time = tmp;
  // fftw_dft(filterf, n, filtert);

  tmp = (complex_t *)calloc(n, sizeof(*tmp));
  memcpy(tmp, filter_est.time, w_est * sizeof(*tmp));
  free(filter_est.time);
  filter_est.time = tmp;
  // fftw_dft(filterf_est, n, filtert_est);

  // plot("filter time series", map_abs(Vec(filtert, n)));
  plot("filter fourier series", map_abs(Vec(filter.freq, 4 * n / B_loc)));

  // plot("filterest time series", map_abs(Vec(filtert_est, n)));
  plot("filterest fourier series",
       map_abs(Vec(filter_est.freq, 4 * n / B_loc)));

  complex_t *out = (complex_t *)malloc(n * sizeof(*out));
  fftw_dft(out, n, x);
  for (int i = 0; i < n; i++)
    out[i] /= n;

  // plot("Original time series", map_real(Vec(x, n)));
  plot("Original fourier series", map_abs(Vec(out, n)));
  free(out);
#endif

  std::map<int, complex_t> ans;
  LOG(INFO) << "n=" << n << " k=" << k;
  reset_timer();
  for (int iii = 0; iii < gInnerReps; iii++) {
    ans =
        outer_loop(x, n, filter, filter_est, B_est, B_thresh, B_loc, W_Comb,
                   Comb_loops, loops_thresh, loops_loc, loops_loc + loops_est);
  }
  SFFT_time = get_time() / double(gInnerReps);

  // JW: Check if we want to compute the error
  if (gComputeError) {

    int num_candidates = (int)ans.size();
    std::pair<real_t, int> *candidates =
        (std::pair<real_t, int> *)malloc(num_candidates * sizeof(*candidates));
    complex_t *x_f_Large = (complex_t *)calloc(n, sizeof(*x_f_Large));
    complex_t *ans_Large = (complex_t *)calloc(n, sizeof(*ans_Large));

    int counter = 0;

    real_t ERROR = 0;
    for (__typeof(ans.begin()) it = ans.begin(); it != ans.end(); it++) {
      int key = it->first;
      complex_t value = it->second;
      candidates[counter] = std::make_pair(cabs(value), key);
      counter++;
    }

    // Enter ALL large frequences as zero
    for (int i = 0; i < k; i++) {
      x_f_Large[LARGE_FREQ[i]] = x_f[LARGE_FREQ[i]];
    }

    std::nth_element(candidates, candidates + num_candidates - k,
                     candidates + num_candidates);
    for (int i = 0; i < k; i++) {
      int key = candidates[num_candidates - k + i].second;
      ans_Large[key] = ans[key];
    }

    int large_found = 0;
    int FOUND = 0;
    for (int i = 0; i < k; i++) {
      FOUND += (unsigned int)ans.count(LARGE_FREQ[i]);
      large_found += (ans_Large[(LARGE_FREQ[i])] != 0);
    }

    // Estimate error as the difference of the K largest entries of x_f and ans)
    for (int i = 0; i < n; i++) {
      ERROR += cabs(ans_Large[i] - x_f_Large[i]);
    }

    error_per_entry = ERROR / k;

    // printf("---------------CONSIDER K LARGEST ONLY
    // ---------------------------------\n");

    // printf("K=%d; MISSED (estimation, result) = (%d, %d); ERROR= %lg  (%lg
    // per entry)\n",k, k-FOUND, k-large_found, ERROR, ERROR/k);

    complex_t *xtmp = (complex_t *)malloc(n * sizeof(*xtmp));
    reset_timer();

    fftw_plan p;

    if (FFTW_OPT)
      p = fftw_plan_dft_1d(n, x, xtmp, FFTW_FORWARD, FFTW_MEASURE);
    else
      p = fftw_plan_dft_1d(n, x, xtmp, FFTW_FORWARD, FFTW_ESTIMATE);

    // printf("Time for FFTW plan: %lf\n", get_time());
    reset_timer();
    fftw_execute(p);
    FFTW_time = get_time();

    // printf("Time for FFTW (%d loops): %lf\n", repetitions, get_time());
    fftw_destroy_plan(p);

    free(xtmp);
    free(candidates);
    free(x_f_Large);
    free(ans_Large);
  }
  free(filter.freq);
  free(filter.time);
  free(filter_est.freq);
  free(filter_est.time);
}

static void usage(const char *progname) {
  const char *p = strrchr(progname, '/'); // drop leading directory path
  if (p)
    p++;
  if (strncmp(p, "lt-", 3) == 0) // drop lt- libtool prefix
    p += 3;

  fprintf(stderr, "Usage: %s [options]\n\n", p);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -h                   show this message and exit\n");
  fprintf(stderr, "  -N Generate Run time Vs N [default]\n");
  fprintf(stderr, "  -K Generate Run time Vs K\n");
  fprintf(stderr, "  -S Generate Error Vs SNR\n");
  fprintf(stderr, "  -R Number of Runs to average [default=10]\n");
  fprintf(stderr, "  -O Use FFTW after optimization\n");
  fprintf(stderr, "  -V Verbose\n");
}

int main(int argc, char **argv) {
  int n = 4 * 128 * 8192;
  int k = 100;
  int repetitions = 10;
  double snr = 100;
  double std_noise = 0;
  double Bcst_loc = 1;
  double Bcst_est = 1;
  double Comb_cst = 2;
  int loc_loops = 4;
  int est_loops = 16;
  int threshold_loops = 3;
  int Comb_loops = 1;
  int ch;
  int Graph_type = 1;
  bool FFTW_OPT = false;
  double tolerance_loc = 1.e-8;
  double tolerance_est = 1.e-8;

  TIMING = false;

  // JW: Replace!
  // int gOuterReps,gInnerReps,gComputeError,gGraphType,gVer,filename;
  if (argc != 7) {
    printf("Need 6 arguments\n");
    exit(0);
  }
  gOuterReps = atoi(argv[1]);
  gInnerReps = atoi(argv[2]);
  gComputeError = atoi(argv[3]);
  gGraphType = atoi(argv[4]);
  gVer = atoi(argv[5]);

  std::string fout_filename(argv[6]);
  CHECK(!fout_filename.empty());
  std::ofstream fout(fout_filename.c_str());
  LOG(INFO) << "Writing to " << fout_filename;

  Graph_type = gGraphType;
  WITH_COMB = (gVer == 2);
  repetitions = gOuterReps;

  /*  while ((ch = getopt(argc, argv, "hNKSR:OWV")) != EOF){
      switch (ch){
      case 'N':
            Graph_type =1;
        break;
      case 'K':
            Graph_type =2;
        break;
      case 'S':
            Graph_type =3;
        break;
          case 'R':
        repetitions = atoi(optarg);
        break;
      case 'O':
        FFTW_OPT = true;
            break;
      case 'W':
        WITH_COMB = true;
            break;
          case 'V':
            VERBOSE = true;
            break;
      case 'h':
      default:
        usage(argv[0]);
        exit(1);
      }
    }*/

  int length = 1;
  int N_vec[12] = {8192,   16384,   32768,   65536,   131072,  262144,
                   524288, 1048576, 2097152, 4194304, 8388608, 16777216};
  int K_vec[8] = {50, 100, 200, 500, 1000, 2000, 2500, 4000};
  double SNR_vec[14] = {-20, -10, -7, -3, 0, 3, 7, 10, 20, 30, 40, 50, 60, 120};

  if (Graph_type == 1)
    length = 12;
  else if (Graph_type == 2)
    length = 8;
  else
    length = 14;

  std::pair<double, double> *SFFT_Time =
      (std::pair<double, double> *)malloc(length * sizeof(*SFFT_Time));
  std::pair<double, double> *FFTW_Time =
      (std::pair<double, double> *)malloc(length * sizeof(*FFTW_Time));
  std::pair<double, double> *SFFT_Error =
      (std::pair<double, double> *)malloc(length * sizeof(*SFFT_Error));

  for (int pp = 0; pp < length; pp++) {

    if (Graph_type == 1) {
      n = N_vec[pp];
      k = 50;
      get_expermient_vs_N_parameters(n, WITH_COMB, Bcst_loc, Bcst_est, Comb_cst,
                                     loc_loops, est_loops, threshold_loops,
                                     Comb_loops, tolerance_loc, tolerance_est);
    } else if (Graph_type == 2) {
      n = 4194304;
      k = K_vec[pp];
      get_expermient_vs_K_parameters(k, WITH_COMB, Bcst_loc, Bcst_est, Comb_cst,
                                     loc_loops, est_loops, threshold_loops,
                                     Comb_loops, tolerance_loc, tolerance_est);
    } else {
      n = 4194304;
      k = 50;
      snr = SNR_vec[pp];
      get_expermient_vs_N_parameters(n, WITH_COMB, Bcst_loc, Bcst_est, Comb_cst,
                                     loc_loops, est_loops, threshold_loops,
                                     Comb_loops, tolerance_loc, tolerance_est);
    }

    // SET FILTER PARAMETERS

    assert(ALGORITHM1 || WITH_COMB);

    real_t BB_loc = (unsigned)(Bcst_loc * sqrt((double)n * k / (log2(n))));
    real_t BB_est = (unsigned)(Bcst_est * sqrt((double)n * k / (log2(n))));

    double lobefrac_loc = 0.5 / (BB_loc);
    double lobefrac_est = 0.5 / (BB_est);
    int b_loc = int(1.2 * 1.1 * ((double)n / BB_loc));
    int b_est = int(1.4 * 1.1 * ((double)n / BB_est));

    int B_loc = floor_to_pow2(BB_loc);
    int B_thresh = 2 * k;
    int B_est = floor_to_pow2(BB_est);

    int W_Comb = floor_to_pow2(Comb_cst * n / B_loc);

    srand(17);
    srand48(time(NULL) ^ (getpid() * 171717));

    if (Graph_type != 3)
      printf("Running SFFT and FFTW %d times for N=%d and K=%d\n", repetitions,
             n, k);
    else
      printf("Running SFFT and FFTW %d times for N=%d and K=%d and SNR=%f dB\n",
             repetitions, n, k, snr);

    double avg_sfft_time = 0;
    double avg_fftw_time = 0;
    double avg_sfft_error = 0;
    double it_sfft;
    double it_fftw;
    double it_error;

    for (int rr = 0; rr < repetitions; rr++) {

      complex_t *x = (complex_t *)malloc(n * sizeof(*x));
      complex_t *x_f = (complex_t *)calloc(n, sizeof(*x_f));
      int *LARGE_FREQ = (int *)malloc(k * sizeof(*LARGE_FREQ));

      // Randomized the None Zero Bins and Generate Time Domain Data
      for (int i = 0; i < k; i++) {
        LARGE_FREQ[i] = (unsigned)floor(drand48() * n);
        x_f[LARGE_FREQ[i]] = 1.0;
      }

      fftw_dft(x, n, x_f, 1);

      if (Graph_type == 3) {
        double snr_achieved;
        std_noise = sqrt(k / (2 * pow(10, snr / 10)));
        snr_achieved = AWGN(x, n, std_noise);
        fftw_dft(x_f, n, x);
        for (int i = 0; i < n; i++)
          x_f[i] /= n;
      }

      run_experiment(x, n, lobefrac_loc, tolerance_loc, b_loc, B_loc, B_thresh,
                     loc_loops, threshold_loops, lobefrac_est, tolerance_est,
                     b_est, B_est, est_loops, W_Comb, Comb_loops, 1, FFTW_OPT,
                     LARGE_FREQ, k, x_f, it_sfft, it_fftw, it_error);

      // JW: Print more
      printf("%d,%d,%.5e,%.5e\n", n, k, it_sfft, it_error);

      avg_sfft_time += it_sfft;
      avg_fftw_time += it_fftw;
      avg_sfft_error += it_error;

      free(x);
      free(x_f);
      free(LARGE_FREQ);
    }

    if (Graph_type == 1) {
      fout << n << "," << k << "," << snr << "," << avg_sfft_time / repetitions << "\n" << std::fflush;
      SFFT_Time[pp] = std::make_pair(n, avg_sfft_time / repetitions);
      FFTW_Time[pp] = std::make_pair(n, avg_fftw_time / repetitions);
    } else if (Graph_type == 2) {
      fout << n << "," << k << "," << snr << "," << avg_sfft_time / repetitions << "\n" << std::fflush;
      SFFT_Time[pp] = std::make_pair(k, avg_sfft_time / repetitions);
      FFTW_Time[pp] = std::make_pair(k, avg_fftw_time / repetitions);
    } else {
      fout << n << "," << k << "," << snr << "," << avg_sfft_error / repetitions << "\n" << std::fflush;
      SFFT_Error[pp] = std::make_pair(snr, avg_sfft_error / repetitions);
    }
  }

  // JW: Don't bother to plot
  // PLOT RESULTS
  /*    std::string plot_options;
      std::string plot_titles;

    if(Graph_type ==1)
        plot_options ="Run time Vs Signal Size (K=50) \n set key top left \n set
    ylabel 'Run Time (sec)' \n set xlabel 'Signal Size N'\n set xrange
    [6000:24000000] \n set xtics ('2^13' 8192, '2^14' 16384, '2^15' 32768,
    '2^16' 65536, '2^17' 131072, '2^18' 262144,  '2^19' 524288, '2^20' 1048576,
    '2^21' 2097152,  '2^22' 4194304, '2^23' 8388608, '2^24' 16777216)\nset
    logscale xy\nz='";
    else if(Graph_type ==2)
        plot_options ="Run time Vs Sparsity (N=2^22) \n set key top left \n set
    ylabel 'Run Time (sec)' \n set xlabel 'Signal Sparsity K'\n set xrange
    [40:5000] \nset logscale xy\nz='";
    else
        plot_options ="Error Vs SNR (N=2^22, K=50) \n set ylabel 'Error per
    Large Frequency' \n set xlabel 'SNR'\n set xrange [-25:130] \nset logscale
    y\nz='";

    plot_titles = "SFFT 1.0\nFFTW";
    if(WITH_COMB)
            plot_titles = "SFFT 2.0\nFFTW";

    if(Graph_type == 3)
            plot(plot_options, plot_titles, Vec(SFFT_Error,length));
    else
            plot(plot_options, plot_titles , Vec(SFFT_Time,length),
    Vec(FFTW_Time,length));
  */

  free(SFFT_Time);
  free(FFTW_Time);

  return 0;
}
