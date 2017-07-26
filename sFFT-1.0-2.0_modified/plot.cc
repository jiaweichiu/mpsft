#include "plot.h"
#include "fftw.h"
#include <cassert>
#include <sstream>

std::vector<real_t> map_abs(std::vector<complex_t> x) {
  std::vector<real_t> y(x.size());
  for (unsigned int i = 0; i < x.size(); i++)
    y[i] = cabs(x[i]);
  return y;
}

std::vector<real_t> map_real(std::vector<complex_t> x) {
  std::vector<real_t> y(x.size());
  for (unsigned int i = 0; i < x.size(); i++)
    y[i] = creal(x[i]);
  return y;
}

std::vector<real_t> map_imag(std::vector<complex_t> x) {
  std::vector<real_t> y(x.size());
  for (unsigned int i = 0; i < x.size(); i++)
    y[i] = cimag(x[i]);
  return y;
}

std::vector<std::string> &split(const std::string &s, char delim,
                                std::vector<std::string> &elems) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}
std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  return split(s, delim, elems);
}

void gnuplot_output(FILE *f, std::string preamble,
                    std::vector<std::vector<std::pair<real_t, real_t> > > plots,
                    std::string titles) {
  unsigned int nump = (unsigned int)(plots.size());
  assert(nump);
  fprintf(f, "%s\n", DEFAULT_PREAMBLE);
  fprintf(f, "%s\n", preamble.c_str());
  fprintf(f, "plot");
  std::vector<std::string> titles_v = split(titles, '\n');
  for (unsigned int j = 0; j < nump; j++) {
    fprintf(f, "%c \"-\"", j ? ',' : ' ');
    if (j < titles_v.size())
      fprintf(f, "title \"%s\"", titles_v[j].c_str());
  }
  fprintf(f, "\n");
  for (unsigned int j = 0; j < nump; j++) {
    for (unsigned int i = 0; i < plots[j].size(); i++) {
      fprintf(f, "%lg %lg", plots[j][i].first, plots[j][i].second);
      fprintf(f, "\n");
    }
    fprintf(f, "e\n");
  }
}

FILE *spawn_gnuplot() {
  int pipefd[2];
  int r = pipe(pipefd);
  assert(!r);
  int cpid = fork();
  if (cpid == 0) {    // child
    close(pipefd[1]); // write end
    dup2(pipefd[0], 0);
    execlp("gnuplot", "gnuplot", "-persist", NULL);
    // execlp("cat", "cat", NULL);
    assert(false);
  }
  close(pipefd[0]); // read end
  FILE *f = fdopen(pipefd[1], "w");
  return f;
}

int plotn(std::string preamble,
          std::vector<std::vector<std::pair<real_t, real_t> > > plots,
          std::string titles) {
  FILE *f = spawn_gnuplot();
  gnuplot_output(f, preamble, plots, titles);
  fclose(f);
  wait(NULL);
  return 0;
}

void plot_fft(std::string preamble, std::vector<complex_t> x, int real) {
  int n = int(x.size());
  std::vector<complex_t> y(n);
  fftw_dft(&y[0], n, &x[0]);
  if (real)
    plot(preamble, map_real(y));
  else
    plot(preamble, map_abs(y));
}

int plotn(std::string preamble, std::vector<std::vector<real_t> > plots,
          std::string titles) {
  std::vector<std::vector<std::pair<real_t, real_t> > > ans;
  for (unsigned int i = 0; i < plots.size(); i++) {
    std::vector<std::pair<real_t, real_t> > v;
    for (unsigned int j = 0; j < plots[i].size(); j++)
      v.push_back(std::make_pair(j, plots[i][j]));
    ans.push_back(v);
  }
  return plotn(preamble, ans, titles);
}
