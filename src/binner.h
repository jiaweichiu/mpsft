#pragma once

#include "base.h"
#include "window.h"

namespace mps {

constexpr int kBinnerSimple = 1;
constexpr int kBinnerFast = 2;

class Binner {
public:
  Binner(const Window &win, Int bits);

  // tau = q +/- (1 << b) where 0 <= b < bits.
  // Num of rows of "out" is 2*bits+1.
  virtual void BinInTime(const CplexArray &x, const Transform &tf, Int q,
                         CplexMatrix *out) = 0;
  virtual void BinInFreq(const ModeMap &mm, const Transform &tf, Int q,
                         CplexMatrix *out) = 0;

  static Binner *Create(int binner_type, const Window &win, Int bits);

protected:
  const Window &win_;
  Int bits_;
  unique_ptr<FFTPlan> plan_;
  unique_ptr<CplexArray> scratch_;
};

class BinnerSimple : public Binner {
public:
  BinnerSimple(const Window &win, Int bits);

  void BinInTime(const CplexArray &x, const Transform &tf, Int q,
                 CplexMatrix *out) override;
  void BinInFreq(const ModeMap &mm, const Transform &tf, Int q,
                 CplexMatrix *out) override;
};

// Make use of symmetry in tau's to half number of sinusoids.
class BinnerFast : public Binner {
public:
  BinnerFast(const Window &win, Int bits);

  void BinInTime(const CplexArray &x, const Transform &tf, Int q,
                 CplexMatrix *out) override;
  void BinInFreq(const ModeMap &mm, const Transform &tf, Int q,
                 CplexMatrix *out) override;

protected:
  unique_ptr<CplexArray> scratch2_;
};

} // namespace mps