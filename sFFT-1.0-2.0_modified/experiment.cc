/*


 */

#include "computefourier.h"
#include "fft.h"
#include "filters.h"
#include "plot.h"
#include "timer.h"
#include "utils.h"
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

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
                    int *LARGE_FREQ, int k, complex_t *x_f) {
  // b = 200;
  // std = 30;
  printf("sFFT filter parameters for: n=%d, k=%d.\n", n, k);
  printf("*********************************************************************"
         "*********\n");
  if (WITH_COMB)
    printf(" Comb Filter: loops: %d mod: %d/%d\n", Comb_loops, B_thresh,
           W_Comb);
  else
    printf(" Comb Filter: none\n");
  if (ALGORITHM1)
    printf(" Location Filter: (numlobes=%.1lf, tol=%lg, b=%d) B: %d/%d loops: "
           "%d/%d\n",
           0.5 / lobefrac_loc, tolerance_loc, b_loc, B_thresh, B_loc,
           loops_thresh, loops_loc);
  else
    printf(" Location Filter: none\n");
  printf(
      " Estimation Filter: (numlobes=%.1lf, tol=%lg, b=%d) B: %d loops: %d\n",
      0.5 / lobefrac_est, tolerance_est, b_est, B_est, loops_est);
  printf("\n");

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
  printf(" Window size: Location Filter : %d; Estimation Filter : %d;\n", w_loc,
         w_est);

  real_t filter_noise = 0, filter_noise_est = 0;
  for (int i = 0; i < 10; i++) {
    filter_noise =
        std::max(filter_noise, std::max(cabs(filter.freq[n / 2 + i]),
                                        cabs(filter.freq[n / 2 - i])));
    filter_noise_est =
        std::max(filter_noise_est, std::max(cabs(filter_est.freq[n / 2 + i]),
                                            cabs(filter_est.freq[n / 2 - i])));
  }
  printf(" Noise in filter: Location Filter : %lg; Estimation Filter %lg\n",
         filter_noise, filter_noise_est);
  printf("*********************************************************************"
         "*********\n\n");

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

  printf("sFFT Results\n");
  printf("*********************************************************************"
         "*********\n");

  std::map<int, complex_t> ans;

  for (int i = 0; i < repetitions; i++) {
    reset_timer();
    ans =
        outer_loop(x, n, filter, filter_est, B_est, B_thresh, B_loc, W_Comb,
                   Comb_loops, loops_thresh, loops_loc, loops_loc + loops_est);
  }

  int num_candidates = (int)ans.size();
  std::pair<real_t, int> *candidates =
      (std::pair<real_t, int> *)malloc(num_candidates * sizeof(*candidates));
  complex_t *x_f_Large = (complex_t *)calloc(n, sizeof(*x_f_Large));
  complex_t *ans_Large = (complex_t *)calloc(n, sizeof(*ans_Large));

  int counter = 0;

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
  real_t ERROR = 0;
  for (int i = 0; i < n; i++) {
    ERROR += cabs(ans_Large[i] - x_f_Large[i]);
  }

  // printf("---------------CONSIDER K LARGEST ONLY
  // ---------------------------------\n");
  printf("ERROR:\n");
  printf("K=%d; MISSED (estimation, result) = (%d, %d); L1 ERROR= %lg  (%lg "
         "per large frequency)\n",
         k, k - FOUND, k - large_found, ERROR, ERROR / k);
  printf("*********************************************************************"
         "*********\n\n");

  printf("FFTW Results\n");
  printf("*********************************************************************"
         "*********\n");
  complex_t *xtmp = (complex_t *)malloc(n * sizeof(*xtmp));
  reset_timer();

  fftw_plan p;

  if (FFTW_OPT)
    p = fftw_plan_dft_1d(n, x, xtmp, FFTW_FORWARD, FFTW_MEASURE);
  else
    p = fftw_plan_dft_1d(n, x, xtmp, FFTW_FORWARD, FFTW_ESTIMATE);

  printf("Time to create FFTW plan: %lf\n", get_time());
  reset_timer();
  for (int i = 0; i < repetitions; i++) {
    fftw_execute(p);
  }
  printf("Time to run FFTW : %lf\n", get_time());
  fftw_destroy_plan(p);
  printf("*********************************************************************"
         "*********\n\n");

  free(xtmp);
  free(candidates);
  free(x_f_Large);
  free(ans_Large);
  free(filter.freq);
  free(filter.time);
  free(filter_est.freq);
  free(filter_est.time);
}

double evaluate_runtime(int n, double lobefrac, double tolerance,
                        double lobefrac2, double tolerance2, int num, int B,
                        int B2, int location_loops, int est_loops,
                        int loop_threshold, int W_Comb, int Comb_loops) {
  int w = int((1 / M_PI) * (1 / lobefrac) * acosh(1. / tolerance));
  int w2 = int((1 / M_PI) * (1 / lobefrac2) * acosh(1. / tolerance2));

  if (WITH_COMB)
    if (num >= W_Comb)
      return -1;
  if (ALGORITHM1) {
    if (num >= B)
      return -1;
    if (loop_threshold > location_loops)
      return -1;
  }
  if (w > n || w2 > n)
    return -1;

  int loops = location_loops + est_loops;
  double m_frac = WITH_COMB ? num * 1. / W_Comb : 1;
  double projected_hits = ALGORITHM1
                              ? binomial_cdf(num * (1. / B - 1. / n),
                                             location_loops, loop_threshold) *
                                        n * m_frac +
                                    num / 2
                              : n * m_frac;
  // XXX B2 for some, B for some
  int k_est = num / 2;
  double projected_noise_on_k =
      2 * binomial_cdf(k_est * (1. / B2 - 1. / n) / 2, loops, (loops + 1) / 2) *
      k_est;
  double projected_error_rate =
      2 * binomial_cdf(k_est * (1. / B2 - 1. / n) / 4, loops, (loops + 1) / 2) *
          n * m_frac +
      projected_noise_on_k;
  // double projected_error_rate = binomial_cdf((num/2) * (1. / B2 - 1./n),
  // est_loops, (est_loops+1)/2) * (projected_hits - num/2);
  printf("Projected error rate: %lg (%lg per large frequency)\n",
         projected_error_rate, projected_noise_on_k);

  double pages_to_set = num * (n / B) * m_frac * location_loops * 1024;
  bool will_array_memset = pages_to_set > n;

  double const_scorearray = n < (1 << 21) ? 0.3 : 1.8;
  double const_permfilt = 38.0;
  double const_Combtime = 90.0;
  double const_estimationtime = WITH_COMB ? 140 : 150;
  double const_grouping = 23;
  double const_group_sort = (WITH_COMB ? 30 : 0);
  double const_bplusctime = 41;

  double time_scorearray = will_array_memset ? const_scorearray * n : 0;
  double time_Comb = const_Combtime * W_Comb * Comb_loops;
  double time_permfilt =
      const_permfilt * (w * 1. * location_loops + w2 * 1. * est_loops);
  double time_grouping =
      (location_loops * (const_grouping * num * (n / B) * m_frac +
                         const_group_sort * num * log(num)) +
       const_scorearray * (!will_array_memset) * pages_to_set);
  double time_estimation = const_estimationtime * projected_hits * loops;
  double time_bplusc = const_bplusctime * (location_loops * B + est_loops * B2);
  double time_total = time_scorearray + time_Comb + time_permfilt +
                      time_grouping + time_estimation + time_bplusc;

  return time_total;
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
  fprintf(stderr, "  -N Number of Frequencies [default= 2^22]\n");
  fprintf(stderr, "  -K NUmber of Energetic Frequencies [default=100]\n");
  fprintf(stderr, "  -R Number of Runs [default=1]\n");
  fprintf(stderr, "  -B Constant for number of Buckets [default=2]\n");
  fprintf(stderr,
          "  -E Constant for number of Buckets for Estimation [default=0.2]\n");
  fprintf(stderr, "  -L Number of estimation loops [default=12]\n");
  fprintf(stderr, "  -l Number of location loops [default=3]\n");
  fprintf(stderr, "  -r Number of required location loops [default=2]\n");
  fprintf(stderr, "  -M Use SFFT2.0 [default=SFFT1.0] and set Constant for "
                  "number of Comb Buckets [default=12]\n");
  fprintf(stderr, "  -m number of times we run Comb [default=1]\n");
  fprintf(stderr, "  -S SNR as ratio [default=inf]\n");
  fprintf(stderr, "  -A USE Alternative algorithm\n");
  fprintf(stderr, "  -O Use FFTW after optimization\n");
  fprintf(stderr, "  -t Tolerance for the Location Filters [default=1.e-6]\n");
  fprintf(stderr,
          "  -e Tolerance for the Estimation Filters [default=1.e-8]\n");
  fprintf(
      stderr,
      "  -s Calculate the expected runtime and error of sFFT and return.\n");
  fprintf(stderr,
          "  -v Verbose : prints detailed timing of each step of the code.\n");
}

int main(int argc, char **argv) {
  int n = 4 * 128 * 8192;
  int k = 100;
  //  int w = 3*(int)(sqrt(k*n));
  int repetitions = 1;
  double Bcst_loc = 2;
  double Bcst_est = 0.2;
  double Comb_cst = 16;
  int loc_loops = 3;
  int est_loops = 12;
  int threshold_loops = 2;
  int Comb_loops = 1;
  int ch;
  int simulate = 0;
  double snr = 1000000000;
  double std_noise = 0;
  bool FFTW_OPT = false;
  double tolerance_loc = 1.e-6;
  double tolerance_est = 1.e-8;

  while ((ch = getopt(argc, argv, "shvN:K:R:B:E:L:l:r:M:m:S:AOt:e:")) != EOF) {
    switch (ch) {

    case 'N':
      n = atoi(optarg);
      break;

    case 'K':
      k = atoi(optarg);
      break;

    case 'R':
      repetitions = atoi(optarg);
      break;
    case 'B':
      Bcst_loc = strtod(optarg, 0);
      break;
    case 'E':
      Bcst_est = strtod(optarg, 0);
      break;
    case 'L':
      est_loops = atoi(optarg);
      break;
    case 'l':
      loc_loops = atoi(optarg);
      if (!loc_loops)
        ALGORITHM1 = false;
      break;
    case 'r':
      threshold_loops = atoi(optarg);
      break;
    case 'M':
      Comb_cst = strtod(optarg, 0);
      WITH_COMB = true;
      break;
    case 'm':
      Comb_loops = atoi(optarg);
      break;
    case 'S':
      snr = strtod(optarg, 0);
      std_noise = sqrt(k / (2 * snr));
      break;
    case 'A':
      ALGORITHM1 = false;
      loc_loops = 0;
      break;
    case 'O':
      FFTW_OPT = true;
      break;
    case 'v':
      VERBOSE = true;
      break;
    case 't':
      tolerance_loc = strtod(optarg, 0);
      break;
    case 'e':
      tolerance_est = strtod(optarg, 0);
      break;
    case 's':
      simulate = 1;
      break;
    case 'h':
    default:
      usage(argv[0]);
      exit(1);
    }
  }
  n = floor_to_pow2(n);
  assert(ALGORITHM1 || WITH_COMB);

  real_t BB_loc = (unsigned)(Bcst_loc * sqrt((double)n * k / (log2(n))));
  real_t BB_est = (unsigned)(Bcst_est * sqrt((double)n * k / (log2(n))));

  double lobefrac_loc = 0.5 / (BB_loc);
  double lobefrac_est = 0.5 / (BB_est);

  int b_loc = int(1.2 * 1.1 * ((double)n / BB_loc));
  // b_loc = 1;
  // real_t BB2 = (unsigned) (Bcst2*sqrt((double)n*k/(log2(n))));
  int b_est = int(1.4 * 1.1 * ((double)n / BB_est));

  int B_loc = floor_to_pow2(BB_loc);
  int B_thresh = 2 * k;
  int B_est = floor_to_pow2(BB_est);

  int W_Comb = floor_to_pow2(Comb_cst * n / B_loc);

  printf("\n\nRUNNING EXPERIMENT: n=%d, k=%d.\n", n, k);

  printf("\n\nSimulation:\n");
  printf("*********************************************************************"
         "*********\n");
  printf("Expected running time: %lg\n",
         evaluate_runtime(n, lobefrac_loc, tolerance_loc, lobefrac_est,
                          tolerance_est, B_thresh, B_loc, B_est, loc_loops,
                          est_loops, threshold_loops, W_Comb, Comb_loops) *
             1e-9);
  printf("*********************************************************************"
         "*********\n");
  if (simulate)
    return 0;

  complex_t *x = (complex_t *)malloc(n * sizeof(*x));

  srand(17);
  srand48(time(NULL) ^ (getpid() * 171717));

  // Randomized the None Zero Bins

  complex_t *x_f = (complex_t *)calloc(n, sizeof(*x_f));

  int *LARGE_FREQ = (int *)malloc(k * sizeof(*LARGE_FREQ));

  // printf("LARGE BINS:");

  //  for(int i = 0; i < repetitions; i++){

  for (int i = 0; i < k; i++) {
    LARGE_FREQ[i] = (unsigned)floor(drand48() * n);
    x_f[LARGE_FREQ[i]] = 1.0; // Will ADD Random Phase and Amplitude Later.
    //      printf("%d, ",LARGE_FREQ[i]);
  }

  printf("\n");

  fftw_dft(x, n, x_f, 1);

  // ADDED NOISE
  double snr_achieved;
  snr_achieved = AWGN(x, n, std_noise);
  if (std_noise != 0)
    printf("SNR = %g / %.2f dB \n\n", snr_achieved, 10 * log10(snr_achieved));

  fftw_dft(x_f, n, x);
  for (int i = 0; i < n; i++)
    x_f[i] /= n;

  run_experiment(x, n, lobefrac_loc, tolerance_loc, b_loc, B_loc, B_thresh,
                 loc_loops, threshold_loops, lobefrac_est, tolerance_est, b_est,
                 B_est, est_loops, W_Comb, Comb_loops, repetitions, FFTW_OPT,
                 LARGE_FREQ, k, x_f);

  free(x);
  free(x_f);
  free(LARGE_FREQ);
  return 0;
}
