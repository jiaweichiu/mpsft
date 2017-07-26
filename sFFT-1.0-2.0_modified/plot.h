#ifndef PLOT_H
#define PLOT_H

#include "fft.h"
#include <stdarg.h>
#include <string>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

#define DEFAULT_PREAMBLE "set style data linespoints"

std::vector<real_t> map_abs(std::vector<complex_t> x);
std::vector<real_t> map_real(std::vector<complex_t> x);
std::vector<real_t> map_imag(std::vector<complex_t> x);

int plotn(std::string preamble,
          std::vector<std::vector<std::pair<real_t, real_t> > > plots,
          std::string titles = "");

int plotn(std::string preamble, std::vector<std::vector<real_t> > plots,
          std::string titles = "");

void plot_fft(std::string preamble, std::vector<complex_t> x, int real = 0);

template <typename T>
inline int
plot(std::string title, std::string titles, std::vector<T> x,
     std::vector<T> y = std::vector<T>(), std::vector<T> z = std::vector<T>(),
     std::vector<T> w = std::vector<T>(), std::vector<T> a = std::vector<T>(),
     std::vector<T> b = std::vector<T>()) {
  std::vector<std::vector<T> > vals;
  vals.push_back(x);
  if (y.size())
    vals.push_back(y);
  if (z.size())
    vals.push_back(z);
  if (w.size())
    vals.push_back(w);
  if (a.size())
    vals.push_back(a);
  if (b.size())
    vals.push_back(b);
  return plotn("set title '" + title + "'", vals, titles);
}

template <typename T>
inline int
plot(std::string title, std::vector<T> x, std::vector<T> y = std::vector<T>(),
     std::vector<T> z = std::vector<T>(), std::vector<T> w = std::vector<T>(),
     std::vector<T> a = std::vector<T>(), std::vector<T> b = std::vector<T>()) {
  return plot(title, "", x, y, z, w, a, b);
}

#endif
