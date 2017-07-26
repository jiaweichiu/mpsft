#include "computefourier.h"
#include "filters.h"

#include "plot.h"
#include "timer.h"
#include "utils.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>

bool WITH_COMB = false;
bool ALGORITHM1 = true;
bool VERBOSE = false;
bool TIMING = true;

#define vprintf(...)                                                           \
  if (VERBOSE) {                                                               \
    printf(__VA_ARGS__);                                                       \
  }

inline const int timesmod(const int &x, const int &a, const int &n) {
  return int((((long long int)x) * a) % n);
}

void debug_inner_loop(const complex_t *origx, int n, const Filter &filter,
                      int num, int B, int a, int ai, int b, int *J,
                      real_t *samples);

int Comb_Filt(const complex_t *origx, int n, int num, int W_Comb,
              int *Comb_Approved) {

  if (n % W_Comb) {
    fprintf(
        stderr,
        "Warning: W_Comb is not divisible by N, which algorithm expects.\n");
    assert(n % W_Comb == 0);
  }

  complex_t *x_sampt = (complex_t *)malloc(W_Comb * sizeof(*x_sampt));

  int sigma = n / W_Comb;
  int offset = (unsigned)floor(drand48() * sigma);

  for (int i = 0; i < W_Comb; i++) {
    x_sampt[i] = origx[offset + i * sigma];
  }

  fftw_dft(x_sampt, W_Comb, x_sampt);
  real_t *samples = (real_t *)malloc(W_Comb * sizeof(*samples));
  for (int i = 0; i < W_Comb; i++)
    samples[i] = cabs2(x_sampt[i]);

  find_largest_indices(Comb_Approved, num, samples, W_Comb);

  free(x_sampt);
  free(samples);
  return 0;
}

/*
  Inner loop of the algorithm, part one.

  n-dimensional origx
  permute the fourier spectrum, take the first w coordinates.
  dot with the filter
  B-dimensional FFT
  return the top num samples.

  Output to
 */
int inner_loop_locate(const complex_t *origx, int n, const Filter &filter,
                      int num, int B, int a, int ai, int b, complex_t *x_samp,
                      int *J, double &PF_T, double &BC_T) {

  if (n % B)
    fprintf(stderr,
            "Warning: n is not divisible by B, which algorithm expects.\n");

  double DDD = get_time(), DDD2 = get_time();
  complex_t *x_sampt = (complex_t *)malloc(n * sizeof(*x_sampt));

  memset(x_sampt, 0, B * sizeof(x_sampt[0]));

  // Permute, dot, collate all in one loop.
  int index = b;
  for (int i = 0; i < filter.sizet; i++) {
    x_sampt[i % B] += origx[index] * filter.time[i];
    index = (index + ai) % n;
  }

  if (TIMING) {
    PF_T = get_time() - DDD2;
    vprintf("Step 1.A (PERM + FILTER):------------------------ %lf\n",
            get_time() - DDD2);
    DDD = get_time();
    DDD2 = DDD;
  }

  fftw_dft(x_samp, B, x_sampt);
  free(x_sampt);

  if (TIMING) {
    vprintf("Step 1.B (FFTW)---------: %lf\n", get_time() - DDD);
    DDD = get_time();
  }

  real_t *samples = (real_t *)malloc(B * sizeof(*samples));
  for (int i = 0; i < B; i++)
    samples[i] = cabs2(x_samp[i]);

  find_largest_indices(J, num, samples, B);

  if (TIMING) {
    BC_T = get_time() - DDD2;
    vprintf("Step 1.C (LARGEST BUCKS): %lf\n", get_time() - DDD);
    DDD = get_time();
    DDD2 = DDD;
  }

#ifdef DEBUG
  debug_inner_loop(origx, n, filter, num, B, a, ai, b, J, samples);
#endif

  free(samples);

  return 0;
}

/*
 Find indices that map to J , i.e., lie within n/(2B) of (J * n/B) after
 permutation.

 For each such i, increment score[i] and append to hits if score[i] reaches
 loop_threshold.
*/
int inner_loop_filter_regular(int *J, int n, int num, int B, int a, int ai,
                              int b, int loop_threshold, int *score, int *hits,
                              int &hits_found, double &G_T) {
  double DDD = get_time();

  // Given the set of large samples, find the locations in [n] that map there
  // and output them

  for (int i = 0; i < num; i++) {
    int low, high;
    low = (int(ceil((J[i] - 0.5) * n / B)) + n) % n;
    high = (int(ceil((J[i] + 0.5) * n / B)) + n) % n;
    int loc = timesmod(low, a, n);
    for (int j = low; j != high; j = (j + 1) % n) {
      score[loc]++;
      if (score[loc] == loop_threshold)
        hits[hits_found++] = loc;
      loc = (loc + a) % n;
    }
  }

  if (TIMING) {
    G_T = get_time() - DDD;
    vprintf("Step 1.D (GROUPING):----------------------------- %lf\n\n",
            get_time() - DDD);
    vprintf("##################################################################"
            "###\n\n");
  }

  return 0;
}

/*
 Find indices that (1) map to J under the permutation and (2) lie in
 Comb_Approved mod W_Comb.

 For each such i, increment hits[i] and append to hits_found if hits[i] reaches
 loop_threshold.

  */
int inner_loop_filter_Comb(int *J, int n, int num, int B, int a, int ai, int b,
                           int loop_threshold, int *score, int *hits,
                           int &hits_found, double &G_T, int *Comb_Approved,
                           int num_Comb, int W_Comb) {
  double DDD = get_time();

  std::pair<int, int> *permuted_approved = (__typeof(permuted_approved))malloc(
      num_Comb * sizeof(*permuted_approved));
  for (int m = 0; m < num_Comb; m++) {
    int prev = timesmod(Comb_Approved[m], ai, W_Comb);
    permuted_approved[m] = std::make_pair(prev, timesmod(prev, a, n));
  }
  std::sort(permuted_approved, permuted_approved + num_Comb);

  // compute intersection of permuted_approved and indices close to J * n/B,
  // then invert to get true locations.

  for (int i = 0; i < num; i++) {
    int low, high;
    low = (int(ceil((J[i] - 0.5) * n / B)) + n) % n;
    high = (int(ceil((J[i] + 0.5) * n / B)) + n) % n;
    int index =
        int(std::upper_bound(permuted_approved, permuted_approved + num_Comb,
                             std::make_pair(low % W_Comb, -1)) -
            permuted_approved);
    int location = low - (low % W_Comb);
    int locinv = timesmod(location, a, n);
    for (int j = index;; j++) {
      if (j == num_Comb) {
        j -= num_Comb;
        location = (location + W_Comb) % n;
        locinv = timesmod(location, a, n);
      }
      int approved_loc = location + permuted_approved[j].first;
      if ((low < high && (approved_loc >= high || approved_loc < low)) ||
          (low > high && (approved_loc >= high && approved_loc < low)))
        break;
      int loc = (locinv + permuted_approved[j].second) % n;
      score[loc]++;
      if (score[loc] == loop_threshold)
        hits[hits_found++] = loc;

      //  printf("{%d,%d}:::(%d-%d)---(%d-%d)---(%d==%d)--:%d--->%d | (%d==%d)||
      //  %d ||| ai=%d, b=%d,
      //  a=%d\n",i,j,low,high,low%W_Comb,high%W_Comb,approved_loc%W_Comb,Comb_Approved[j],approved_loc,loc,loc%W_Comb,mod_approved,score[loc],ai,b,a);
    }
  }
  free(permuted_approved);

  if (TIMING) {
    G_T = get_time() - DDD;
    vprintf("Step 1.D (GROUPING):----------------------------- %lf\n\n",
            get_time() - DDD);
    vprintf("##################################################################"
            "###\n\n");
  }

  return 0;
}

/*
  hits contains the indices that we want to estimate.

  x_samp contains a B-dimensional array for each of the `loops`
  iterations of the outer loop.  Every coordinate i of x "hashes to" a
  corresponding coordinate (permute[j] * i) mod B of x_samp[j], which
  gives an estimate of x[i].

  We estimate each coordinate as the median (independently in real and
  imaginary axes) of its images in the rows of x_samp.
 */
std::map<int, complex_t> estimate_values(const int *hits, const int &hits_found,
                                         complex_t **x_samp, const int &loops,
                                         int n, const int *permute, const int B,
                                         const int B2, const Filter &filter,
                                         const Filter &filter_Est,
                                         int location_loops) {
  std::map<int, complex_t> ans;
  real_t *values[2];
  for (int a = 0; a < 2; a++)
    values[a] = (real_t *)malloc(loops * sizeof(*values[a]));

  for (int i = 0; i < hits_found; i++) {
    int position = 0;

    for (int j = 0; j < loops; j++) {
      int cur_B = (j < location_loops) ? B : B2;
      const Filter &cur_filter = (j < location_loops) ? filter : filter_Est;
      int permuted_index = timesmod(permute[j], hits[i], n);
      int hashed_to = permuted_index / (n / cur_B);
      int dist = permuted_index % (n / cur_B);
      if (dist > (n / cur_B) / 2) {
        hashed_to = (hashed_to + 1) % cur_B;
        dist -= n / cur_B;
      }
      dist = (n - dist) % n;
      complex_t filter_value = cur_filter.freq[dist]; // * cexp(2*M_PI * I *
                                                      // timesmod(permuteb[j],
                                                      // hits[i], n) / n);
      values[0][position] = creal(x_samp[j][hashed_to] / filter_value);
      values[1][position] = cimag(x_samp[j][hashed_to] / filter_value);
      position++;
      // printf("MOO %d %lf %lf: %lf %d %lf %lf+%lfj\n", hits[i], permuted_index
      // * 1./n, hashed_to * 1./B, hashed_to * (n * 1. /B), dist,
      // cabs(filter_value), values[0][position-1], values[1][position-1]);
    }

    int location = (loops - 1) / 2;

    for (int a = 0; a < 2; a++)
      std::nth_element(values[a], values[a] + location, values[a] + position);
    real_t realv = values[0][location];
    real_t imagv = values[1][location];
    ans[hits[i]] = realv + I * imagv;
  }
  for (int a = 0; a < 2; a++)
    free(values[a]);
  return ans;
}

/*
  Outer loop of the algorithm.

  If we are performing the Comb heuristic, first we do so.

  Then, `loops` times:
    choose a random permutation
    run inner_loop_locate
    if in the first location_loops loops, also run inner_loop_filter

  at the end, `hits` contains the coordinates that appear at least
  loop_threshold of location_loops times.  We estimate the values at
  these coordinates as the median of the images x_samp[loops].

  Returns a map from coordinates to estimates.
 */
std::map<int, complex_t>
outer_loop(const complex_t *origx, int n, const Filter &filter,
           const Filter &filter_Est, int B2, int num, int B, int W_Comb,
           int Comb_loops, int loop_threshold, int location_loops, int loops) {
  int *permute = (int *)malloc(loops * sizeof(*permute));
  int *permuteb = (int *)malloc(loops * sizeof(*permuteb));

  complex_t *x_samp[loops];
  for (int i = 0; i < loops; i++) {
    if (i < location_loops)
      x_samp[i] = (complex_t *)calloc(B, sizeof(*x_samp[i]));
    else
      x_samp[i] = (complex_t *)calloc(B2, sizeof(*x_samp[i]));
  }
  int hits_found = 0;

  // Variables used for timing
  double SCORE_T = 0;
  double PF_T = 0;
  double G_T = 0;
  double BC_T = 0;
  double PF_ALL = 0;
  double G_ALL = 0;
  double BC_ALL = 0;

  double DDD = get_time();

  // calloc is faster if few pages are hit, while malloc/memset is
  // faster if most pages are hit.
  double pages_hit =
      num * (n / B) * (WITH_COMB ? num * 1. / W_Comb : 1.) * location_loops;
  int *score;
  if (pages_hit > n / 1024) {
    score = (int *)malloc(n * sizeof(*score));
    memset(score, 0, n * sizeof(*score));
  } else {
    score = (int *)calloc(n, sizeof(*score));
  }

  SCORE_T = get_time() - DDD;
  // printf("Created score array: %lf\n", get_time() - DDD);

  int *hits = (int *)malloc(n * sizeof(*hits));

  double PF_LOC = 0;
  double G_LOC = 0;

  // BEGIN Comb
  DDD = get_time();
  int *Comb_Approved = (int *)malloc(Comb_loops * num * sizeof(*Comb_Approved));
  int num_Comb = num;

  if (WITH_COMB) {
    for (int i = 0; i < Comb_loops; i++)
      Comb_Filt(origx, n, num, W_Comb, Comb_Approved + i * num);
  }

  if (Comb_loops > 1) {
    radix_sort(Comb_Approved, Comb_loops * num);
    int Last = 0;
    for (int i = 1; i < Comb_loops * num; i++) {
      if (Comb_Approved[i] != Comb_Approved[Last])
        Comb_Approved[++Last] = Comb_Approved[i];
    }
    num_Comb = Last + 1;

    vprintf("Comb:%d----->%d\n\n", num * Comb_loops, num_Comb);
  }

  if (!ALGORITHM1) {
    hits_found = num_Comb * (n / W_Comb);
    for (int j = 0; j < n / W_Comb; j++)
      for (int i = 0; i < num_Comb; i++)
        hits[j * num_Comb + i] = j * W_Comb + Comb_Approved[i];
  }

  double Comb_time = get_time() - DDD;
  // END Comb

  // BEGIN INNER LOOPS

  for (int i = 0; i < loops; i++) {
    int a = 0;
    int b = 0; // random() % n;
    while (gcd(a, n) != 1) {
      a = int(random() % n);
    }
    int ai = mod_inverse(a, n);

    permute[i] = ai;
    permuteb[i] = b;

    int perform_location = (i < location_loops);
    assert(ALGORITHM1 || !perform_location);
    Filter cur_filter = perform_location ? filter : filter_Est;
    int cur_B = perform_location ? B : B2;

    int *J = (int *)malloc(num * sizeof(*J));

    inner_loop_locate(origx, n, cur_filter, num, cur_B, a, ai, b, x_samp[i], J,
                      PF_T, BC_T);
    if (perform_location) {
      if (!WITH_COMB) {
        inner_loop_filter_regular(J, n, num, cur_B, a, ai, b, loop_threshold,
                                  score, hits, hits_found, G_T);
      } else {
        inner_loop_filter_Comb(J, n, num, cur_B, a, ai, b, loop_threshold,
                               score, hits, hits_found, G_T, Comb_Approved,
                               num_Comb, W_Comb);
      }
    }
    free(J);

    PF_ALL += PF_T;
    BC_ALL += BC_T;
    if (perform_location) {
      PF_LOC += PF_T;
      G_LOC += G_T;
      G_ALL += G_T;
    }
  }

  vprintf("Number of candidates: %d\n", hits_found);
  // END INNER LOOPS

  // BEGIN ESTIMATION
  DDD = get_time();

  std::map<int, complex_t> ans =
      estimate_values(hits, hits_found, x_samp, loops, n, permute, B, B2,
                      filter, filter_Est, location_loops);

  // END ESTIMATION
  double E_T = get_time() - DDD;

#ifdef DEBUG
  real_t *x = (real_t *)malloc(n * sizeof(*x));
  memset(x, 0, n * sizeof(*x));
  for (std::map<int, complex_t>::iterator it = ans.begin(); it != ans.end();
       it++) {
    x[it->first] = cabs(it->second);
  }

  real_t *xc = (real_t *)malloc(n * sizeof(*xc));
  memset(xc, 0, n * sizeof(*xc));
  for (int i = 0; i < n; i++)
    xc[i] = score[i] * 1. / loops;

  complex_t *xf = (complex_t *)malloc(n * sizeof(*xf));
  fftw_dft(xf, n, origx);
  for (int i = 0; i < n; i++)
    xf[i] /= n;
  plot("outer loop", "counts\ntrue signal\nreconstruction", Vec(xc, n),
       map_abs(Vec(xf, n)), Vec(x, n));

  free(x);
  free(xc);
  free(xf);
#endif // DEBUG

  for (int i = 0; i < loops; i++)
    free(x_samp[i]);
  free(permute);
  free(score);
  free(hits);
  free(Comb_Approved);

  DDD = get_time();

  if (TIMING) {
    printf("Total sFFT time: %lf\n", DDD);
    printf("Time distribution: scoretable  Comb __  perm+filter grouping "
           "estimation  stepB+C    other    total\n");
    printf("                     %lf %lf    %lf %lf   %lf %lf %lf %lf\n",
           SCORE_T, Comb_time, PF_ALL, G_ALL, E_T, BC_ALL,
           DDD - PF_ALL - G_ALL - E_T - Comb_time - BC_ALL - SCORE_T, DDD);
    double tott = (DDD) / 100;
    printf("                        %4.1lf%%    %4.1lf%%       %4.1lf%%    "
           "%4.1lf%%      %4.1lf%%    %4.1lf%%    %4.1lf%%   %5.1lf%%\n",
           SCORE_T / tott, Comb_time / tott, PF_ALL / tott, G_ALL / tott,
           E_T / tott, BC_ALL / tott,
           (DDD - PF_ALL - G_ALL - E_T - Comb_time - BC_ALL - SCORE_T) / tott,
           (DDD) / tott);

    // printf("LOC/EST loops: perm+filter grouping total\n");
    // printf("Location:         %lf %lf %lf\n", PF_LOC, G_LOC, PF_LOC + G_LOC);
    // printf("Estimation:       %lf %lf %lf\n", PF_ALL-PF_LOC, G_ALL-G_LOC,
    // (PF_ALL + G_ALL) - (PF_LOC + G_LOC));
    printf("\n");
  }

  return ans;
}

///////////////////////////////////////////////////////////////////
// Everything beyond this point is for debugging / testing only. //
///////////////////////////////////////////////////////////////////

/*

  predict_performance uses empirically determined constants to predict
  the running time for a given set of parameters.  In principle, we
  could use this to optimize the choice of parameters.


  On the computer for which the constants were tuned, it is reliably
  within 10% of the true running time.  Unfortunately, the constants
  vary substantially by computer.

 */

double predict_performance(int n, const Filter &filter,
                           const Filter &filter_Est, int B2, int num, int B,
                           int W_Comb, int Comb_loops, int loop_threshold,
                           int location_loops, int loops) {
  int est_loops = loops - location_loops;
  double scale = 1e9;
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
  // double const_estimationtime = WITH_COMB ? 210 : 220;
  double const_estimationtime = WITH_COMB ? 140 : 150;
  double const_grouping = 23;
  double const_group_sort = (WITH_COMB ? 30 : 0);
  double const_bplusctime = 41;

  double time_scorearray = will_array_memset ? const_scorearray * n : 0;
  double time_Comb = const_Combtime * W_Comb * Comb_loops;
  double time_permfilt = const_permfilt * (filter.sizet * 1. * location_loops +
                                           filter_Est.sizet * 1. * est_loops);
  double time_grouping =
      (location_loops * (const_grouping * num * (n / B) * m_frac +
                         const_group_sort * num * log(num)) +
       const_scorearray * (!will_array_memset) * pages_to_set);
  double time_estimation = const_estimationtime * projected_hits * loops;
  double time_bplusc = const_bplusctime * (location_loops * B + est_loops * B2);
  double time_total = time_scorearray + time_Comb + time_permfilt +
                      time_grouping + time_estimation + time_bplusc;

  printf("Predictions:       scoretable  Comb perm+filter grouping estimation  "
         "stepB+C             total\n");
  printf("                     %lf %lf    %lf %lf   %lf %lf          %lf\n",
         time_scorearray / scale, time_Comb / scale, time_permfilt / scale,
         time_grouping / scale, time_estimation / scale, time_bplusc / scale,
         time_total / scale);
  //    printf("real/estimate        %7.1lf%% %7.1lf%%    %7.1lf%% %7.1lf%%
  //    %7.1lf%% %7.1lf%%          %7.1lf%%\n",
  //       100*scale *SCORE_T/time_scorearray, 100*scale *Comb_time/time_Comb,
  //       100*scale *PF_ALL/time_permfilt, 100*scale *G_ALL/time_grouping,
  //       100*scale *E_T/time_estimation, 100*scale *BC_ALL/time_bplusc,
  //       100*scale*DDD/time_total);

  printf("Projected hits_found: %d\n", int(projected_hits));

  return time_total;
}

/*

  When called in the inner loop, this function plots some informative graphics
  about what's going on.

 */
void debug_inner_loop(complex_t *origx, int n, const Filter &filter, int num,
                      int B, int a, int ai, int b, int *J, real_t *samples) {
  complex_t *xf = (complex_t *)malloc(n * sizeof(*xf));
  complex_t *pxdotg = (complex_t *)malloc(n * sizeof(*pxdotg));
  complex_t *pxdotgn = (complex_t *)malloc(n * sizeof(*pxdotg));
  complex_t *pxdotgw = (complex_t *)malloc(n * sizeof(*pxdotg));
  memset(xf, 0, n * sizeof(*xf));
  int index = b;
  for (int i = 0; i < n; i++) {
    xf[i] = origx[index];
    index = (index + ai) % n;
  }

  int w = filter.sizet;
  memcpy(pxdotg, xf, w * sizeof(*pxdotg));
  memset(pxdotg + w, 0, (n - w) * sizeof(*pxdotg));
  for (int i = 0; i < w; i++)
    pxdotg[i] *= filter.time[i];
  // plot("time after permutation", map_real(Vec(xf, n)));
  // plot_fft("FFT after permutation", Vec(xf, n));
  printf("Using %dx (%d^-1)\n", a, ai);
  fftw_dft(xf, n, xf);

  fftw_dft(pxdotgn, n, pxdotg);
  fftw_dft(pxdotgw, w, pxdotg);

  real_t xtmp[B];
  memset(xtmp, 0, B * sizeof(*xtmp));
  for (int i = 0; i < num; i++)
    xtmp[J[i]] = sqrt(samples[J[i]]);
  std::vector<std::pair<real_t, real_t> > fft_n, fft_w, fft_sample, fft_large,
      fft_true;
  for (int i = 0; i < n; i++) {
    fft_n.push_back(std::make_pair(i * 1. / n, cabs(pxdotgn[i])));
    fft_true.push_back(std::make_pair(i * 1. / n, cabs(xf[i]) / n));
  }
  for (int i = 0; i < w; i++)
    fft_w.push_back(std::make_pair(i * 1. / w, cabs(pxdotgw[i])));
  for (int i = 0; i < B; i++) {
    fft_sample.push_back(std::make_pair(i * 1. / B, sqrt(samples[i])));
  }
  for (int i = 0; i < num; i++) {
    if (J[i])
      fft_large.push_back(std::make_pair((J[i] - 0.1) * 1. / B, 0));
    fft_large.push_back(std::make_pair(J[i] * 1. / B, sqrt(samples[J[i]])));
    if (J[i] < B - 1)
      fft_large.push_back(std::make_pair((J[i] + 0.1) * 1. / B, 0));
  }
  plot("Inner loop result", "n-dim convolved FFT\nw-dim convolved FFT\nsampled "
                            "convolved FFT\nlargest in sample\ntrue FFT",
       fft_n, fft_w, fft_sample, fft_large, fft_true);
  free(xf);
  free(pxdotg);
  free(pxdotgn);
  free(pxdotgw);
}
