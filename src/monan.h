// Data structures of the C++ MCMC sampler: one chain (edgelist and cache) and
// the interface of effects. Indices are 0-based; matrices are column-major
// n x n vectors, as in R.
#ifndef MONAN_H
#define MONAN_H

// No fused multiply-add: R rounds a + b * c in two steps, and the C++ effects
// must give exactly the same values as the R effects. (A pragma rather than
// -ffp-contract=off in Makevars, which R CMD check reports as non-portable.)
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include <Rcpp.h>
#include <Rmath.h>
#include <memory>
#include <string>
#include <vector>

// R's x^y (arithmetic.c): x * x for y == 2, R_pow() otherwise. R_pow() is
// R's own power function, which does not use pow() on all platforms (on
// 64-bit Windows, it uses the long double powl()).
inline double rpow(double x, double y) {
  if (y == 2.0) return x * x;
  return R_pow(x, y);
}

struct Chain {
  int n;                                  // number of locations
  std::vector<int> sender, receiver;      // per resource; -1 while it is moved
  std::vector<std::string> resourceNames; // resource covariates in the cache
  std::vector<std::vector<double>> resourceValues;  // per covariate, per resource

  // cache (as cache[[dep.var]] in R)
  std::vector<double> valued, netFlows, minNet;
  std::vector<double> twoPaths;           // twoPathNetwork (getTwoPathNetwork in R)
  std::vector<std::vector<double>> resourceNetworks;

  inline int idx(int i, int j) const { return i + j * n; }
  inline double v(int i, int j) const { return valued[idx(i, j)]; }
  inline double tp(int i, int j) const { return twoPaths[idx(i, j)]; }

  // same logic as updateWeightedCache / the in-place .update in R
  void update(int s, int r, int edge, int u) {
    for (size_t k = 0; k < resourceNetworks.size(); ++k) {
      resourceNetworks[k][idx(s, r)] += u * resourceValues[k][edge];
    }
    // the following uses information about the net flow network before the update
    if (u > 0 && netFlows[idx(s, r)] < 0) minNet[idx(s, r)] += 1;
    if (u < 0 && netFlows[idx(s, r)] <= 0 && s != r) minNet[idx(s, r)] -= 1;
    if (minNet[idx(s, r)] < 0) {
      Rcpp::stop("Error in cache update, negative min tie. %d %d %d", s + 1, r + 1, u);
    }
    minNet[idx(r, s)] = minNet[idx(s, r)];
    // two-paths, as twoPathIncrements in R: row s gains u times row r and
    // column r gains u times column s, both without loops
    if (s != r) {
      for (int k = 0; k < n; ++k) {
        if (k != r) twoPaths[idx(s, k)] += u * v(r, k);
        if (k != s) twoPaths[idx(k, r)] += u * v(k, s);
      }
    }
    valued[idx(s, r)] += u;
    netFlows[idx(s, r)] += u;
    netFlows[idx(r, s)] -= u;
  }

  // cells changed by an update of s-r, to restore them exactly
  struct Cells {
    double valued, nf_sr, nf_rs, min_sr, min_rs;
    std::vector<double> resource;
    std::vector<double> tpRow, tpCol;     // row s and column r of twoPaths
  };
  Cells getCells(int s, int r) const {
    Cells c{v(s, r), netFlows[idx(s, r)], netFlows[idx(r, s)],
            minNet[idx(s, r)], minNet[idx(r, s)], {}};
    for (const auto& m : resourceNetworks) c.resource.push_back(m[idx(s, r)]);
    c.tpRow.resize(n);
    c.tpCol.resize(n);
    for (int k = 0; k < n; ++k) {
      c.tpRow[k] = tp(s, k);
      c.tpCol[k] = tp(k, r);
    }
    return c;
  }
  void restoreCells(int s, int r, const Cells& c) {
    for (size_t k = 0; k < resourceNetworks.size(); ++k) {
      resourceNetworks[k][idx(s, r)] = c.resource[k];
    }
    valued[idx(s, r)] = c.valued;
    netFlows[idx(s, r)] = c.nf_sr;
    netFlows[idx(r, s)] = c.nf_rs;
    minNet[idx(s, r)] = c.min_sr;
    minNet[idx(r, s)] = c.min_rs;
    for (int k = 0; k < n; ++k) {
      twoPaths[idx(s, k)] = c.tpRow[k];
      twoPaths[idx(k, r)] = c.tpCol[k];
    }
  }

  int resourceIndex(const std::string& name) const {
    for (size_t k = 0; k < resourceNames.size(); ++k) {
      if (resourceNames[k] == name) return (int)k;
    }
    Rcpp::stop("resource covariate '%s' is not in the cache", name);
  }
};

// An effect: change statistic of adding (update = 1) or removing (update = -1)
// resource `edge` on i -> j, and target contribution of the i -> j cell.
class Effect {
 public:
  virtual ~Effect() {}
  virtual double change(const Chain& ch, int i, int j, int edge, int update) const = 0;
  virtual double target(const Chain& ch, int i, int j) const = 0;
  // the target statistic: the sum of the target contributions over all cells,
  // in the order of getNetworkStatistics() (i fastest); effects can override
  // this with a faster calculation that gives exactly the same value
  virtual double targetTotal(const Chain& ch) const {
    long double total = 0;
    for (int j = 0; j < ch.n; ++j) {
      for (int i = 0; i < ch.n; ++i) total += target(ch, i, j);
    }
    return (double)total;
  }
};

// builds an effect from its specification (see cppEffectSpec in R)
std::unique_ptr<Effect> makeEffect(const Rcpp::List& spec, const Chain& chain);

#endif
