// C++ versions of the effects, translated line by line from the R effect
// functions (R/2*_effectFunctions*.R). Indices are 0-based.
#include "monan.h"

namespace {

// loops (R/2b_effectFunctions_loops.R)
class Loops : public Effect {
 public:
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (i == j) return update;
    return 0;
  }
  double target(const Chain& ch, int i, int j) const override {
    return ch.v(i, j) * (i == j);
  }
};

// marginal and cumulative contributions of the AC effects
inline double gMarAC(double y, double a) {
  if (y > 0) return 1 - rpow(1 - 1 / a, y);
  return 0;
}
inline double gCumAC(double y, double a) {
  double contr = 0;
  if (y > 1) {
    for (int k = 1; k <= (int)(y - 1); ++k) contr = contr + (1 - rpow(1 - 1 / a, k));
  }
  return contr;
}

double alphaOf(const Rcpp::List& spec, const char* effect) {
  double alpha = Rcpp::as<double>(spec["alpha"]);
  if (alpha < 1) Rcpp::stop("alpha parameter in %s function must be 1 or larger", effect);
  return alpha;
}

std::vector<double> attributeOf(const Rcpp::List& spec) {
  Rcpp::NumericVector a = spec["attribute"];
  return std::vector<double>(a.begin(), a.end());
}

bool loopContributionOf(const Rcpp::List& spec) {
  return Rcpp::as<bool>(spec["loopContribution"]);
}

// concentration_AC (R/2c_effectFunctions_concentration.R)
class ConcentrationAC : public Effect {
  double alpha;
 public:
  explicit ConcentrationAC(double a) : alpha(a) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (i == j) return 0;
    double tie = ch.v(i, j);
    if (update < 0) return update * gMarAC(tie + update, alpha);
    return update * gMarAC(tie, alpha);
  }
  double target(const Chain& ch, int i, int j) const override {
    if (i == j) return 0;
    return gCumAC(ch.v(i, j), alpha);
  }
};

// reciprocity_AC (R/2d_effectFunctions_reciprocity.R)
class ReciprocityAC : public Effect {
  double alpha;
 public:
  explicit ReciprocityAC(double a) : alpha(a) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (i == j) return 0;
    double nIncoming = ch.v(j, i);
    double cont = update * gMarAC(nIncoming, alpha);
    double nOutgoing = ch.v(i, j);
    cont = cont + nIncoming * (gMarAC(nOutgoing + update, alpha) - gMarAC(nOutgoing, alpha));
    return cont;
  }
  double target(const Chain& ch, int i, int j) const override {
    if (i == j) return 0;
    return ch.v(i, j) * gMarAC(ch.v(j, i), alpha);
  }
};

// dyadic_covariate (R/2a_effectFunctions_exogenous.R)
class DyadicCovariate : public Effect {
  std::vector<double> attribute;  // n x n, column-major
 public:
  explicit DyadicCovariate(std::vector<double> a) : attribute(std::move(a)) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (i == j) return 0;
    return update * attribute[ch.idx(i, j)];
  }
  double target(const Chain& ch, int i, int j) const override {
    if (i == j) return 0;
    return ch.v(i, j) * attribute[ch.idx(i, j)];
  }
};

// alter_covariate (R/2a_effectFunctions_exogenous.R)
class AlterCovariate : public Effect {
  std::vector<double> attribute;
  bool loopContribution;
 public:
  AlterCovariate(std::vector<double> a, bool loops) : attribute(std::move(a)), loopContribution(loops) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (!loopContribution && i == j) return 0;
    return attribute[j] * update;
  }
  double target(const Chain& ch, int i, int j) const override {
    if (!loopContribution && i == j) return 0;
    return attribute[j] * ch.v(i, j);
  }
};

// resource_covar_to_node_covar (R/2a_effectFunctions_exogenous.R)
class ResourceCovarToNodeCovar : public Effect {
  int resource;
  std::vector<double> attribute;
  bool loopContribution;
 public:
  ResourceCovarToNodeCovar(int r, std::vector<double> a, bool loops)
    : resource(r), attribute(std::move(a)), loopContribution(loops) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (!loopContribution && i == j) return 0;
    return update * ch.resourceValues[resource][edge] * attribute[j];
  }
  double target(const Chain& ch, int i, int j) const override {
    if (!loopContribution && i == j) return 0;
    return attribute[j] * ch.resourceNetworks[resource][ch.idx(i, j)];
  }
};

// loops_resource_covar (R/2b_effectFunctions_loops.R)
class LoopsResourceCovar : public Effect {
  int resource;
 public:
  explicit LoopsResourceCovar(int r) : resource(r) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (i != j) return 0;
    return update * ch.resourceValues[resource][edge];
  }
  double target(const Chain& ch, int i, int j) const override {
    if (i != j) return 0;
    return ch.resourceNetworks[resource][ch.idx(i, j)];
  }
};

// loops_AC (R/2b_effectFunctions_loops.R)
class LoopsAC : public Effect {
  double alpha;
 public:
  explicit LoopsAC(double a) : alpha(a) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (i != j) return 0;
    double tie = ch.v(i, i);
    if (update < 0) return update * gMarAC(tie + update, alpha);
    return update * gMarAC(tie, alpha);
  }
  double target(const Chain& ch, int i, int j) const override {
    if (i == j) return gCumAC(ch.v(i, i), alpha);
    return 0;
  }
};

// R's 0:y: counts up to floor(y) for y >= 0, and down to ceiling(y) otherwise
inline std::vector<int> rSeqFromZero(double y) {
  std::vector<int> k;
  if (y >= 0) {
    for (int m = 0; m <= (int)std::floor(y); ++m) k.push_back(m);
  } else {
    for (int m = 0; m >= (int)std::ceil(y); --m) k.push_back(m);
  }
  return k;
}

// marginal and cumulative contributions of the GW effects, summed term by term
// over 0:y as in R
inline double gMarGW(double y, double a) {
  double contr = 0;
  for (int k : rSeqFromZero(y)) contr = contr + std::exp(-std::log(a) * k);
  return contr - 1;
}
inline double gCumGW(double y, double a) {
  double contr = 0;
  for (int k : rSeqFromZero(y)) contr = contr + (y - k) * std::exp(-std::log(a) * k);
  return contr - y;
}

// loops_GW (R/2b_effectFunctions_loops.R)
class LoopsGW : public Effect {
  double alpha;
 public:
  explicit LoopsGW(double a) : alpha(a) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (i != j) return 0;
    double tie = ch.v(i, i);
    if (update < 0) return update * gMarGW(tie + update, alpha);
    return update * gMarGW(tie, alpha);
  }
  double target(const Chain& ch, int i, int j) const override {
    if (i == j) return gCumGW(ch.v(i, i), alpha);
    return 0;
  }
};

// loops_additional_origin and loops_x_loops_additional_origin
// (R/2b_effectFunctions_loops.R): resources moving to their additional origin
// (a resource covariate holding a location id); with onlyLoops, only if this
// is also their origin
class LoopsAdditionalOrigin : public Effect {
  int resource;
  bool onlyLoops;
  bool counts(const Chain& ch, int e, int i, int j) const {
    return ch.sender[e] == i && ch.receiver[e] == j &&
           ch.resourceValues[resource][e] == j + 1;
  }
 public:
  LoopsAdditionalOrigin(int r, bool loops) : resource(r), onlyLoops(loops) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (j + 1 == ch.resourceValues[resource][edge] && (!onlyLoops || j == i)) return update;
    return 0;
  }
  double target(const Chain& ch, int i, int j) const override {
    if (onlyLoops && i != j) return 0;
    double cont = 0;
    for (size_t e = 0; e < ch.sender.size(); ++e) cont += counts(ch, (int)e, i, j);
    return cont;
  }
  // one pass over the resources instead of one per cell; the same count
  double targetTotal(const Chain& ch) const override {
    double total = 0;
    for (size_t e = 0; e < ch.sender.size(); ++e) {
      int i = ch.sender[e], j = ch.receiver[e];
      if (i < 0 || (onlyLoops && i != j)) continue;
      total += ch.resourceValues[resource][e] == j + 1;
    }
    return total;
  }
};

// loops_node_covar (R/2b_effectFunctions_loops.R)
class LoopsNodeCovar : public Effect {
  std::vector<double> attribute;
 public:
  explicit LoopsNodeCovar(std::vector<double> a) : attribute(std::move(a)) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (i != j) return 0;
    return update * attribute[i];
  }
  double target(const Chain& ch, int i, int j) const override {
    if (i != j) return 0;
    return ch.v(i, j) * attribute[i];
  }
};

// loops_resource_covar_node_covar (R/2b_effectFunctions_loops.R)
class LoopsResourceCovarNodeCovar : public Effect {
  int resource;
  std::vector<double> attribute;
 public:
  LoopsResourceCovarNodeCovar(int r, std::vector<double> a) : resource(r), attribute(std::move(a)) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (i == j) return update * ch.resourceValues[resource][edge] * attribute[j];
    return 0;
  }
  double target(const Chain& ch, int i, int j) const override {
    return (ch.resourceNetworks[resource][ch.idx(i, j)] * (i == j)) * attribute[i];
  }
};

// concentration_GW, concentration_AC_dyad_covar, concentration_GW_dyad_covar,
// concentration_AC_resource_covar_bin, concentration_GW_resource_covar_bin
// (R/2c_effectFunctions_concentration.R): the AC or GW function of the tie
// value, which is the number of resources (weighted by a dyadic covariate) or
// the sum of a resource covariate on the tie
class Concentration : public Effect {
 public:
  enum Weight { NONE, DYAD, RESOURCE };
 private:
  bool gw;
  Weight weight;
  double alpha;
  std::vector<double> dyad;  // DYAD: n x n covariate
  int resource;              // RESOURCE: resource covariate
  double tieValue(const Chain& ch, int i, int j) const {
    if (weight == RESOURCE) return ch.resourceNetworks[resource][ch.idx(i, j)];
    return ch.v(i, j);
  }
  double gMar(double y) const { return gw ? gMarGW(y, alpha) : gMarAC(y, alpha); }
  double gCum(double y) const { return gw ? gCumGW(y, alpha) : gCumAC(y, alpha); }
 public:
  Concentration(bool isGW, Weight w, double a, std::vector<double> d, int r)
    : gw(isGW), weight(w), alpha(a), dyad(std::move(d)), resource(r) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (i == j) return 0;
    if (weight == DYAD && dyad[ch.idx(i, j)] == 0) return 0;
    if (weight == RESOURCE && ch.resourceValues[resource][edge] == 0) return 0;
    double tie = tieValue(ch, i, j);
    double g = update < 0 ? update * gMar(tie + update) : update * gMar(tie);
    if (weight == DYAD) return g * dyad[ch.idx(i, j)];
    return g;
  }
  double target(const Chain& ch, int i, int j) const override {
    if (i == j) return 0;
    double g = gCum(tieValue(ch, i, j));
    if (weight == DYAD) return g * dyad[ch.idx(i, j)];
    return g;
  }
};

// reciprocity_GW, reciprocity_AC_dyad_covar, reciprocity_AC_dyad_covar_bin,
// reciprocity_GW_dyad_covar, reciprocity_GW_dyad_covar_bin
// (R/2d_effectFunctions_reciprocity.R)
class Reciprocity : public Effect {
 public:
  enum Weight { NONE, DYAD, DYAD_BIN };
 private:
  bool gw;
  Weight weight;
  double alpha;
  std::vector<double> dyad;  // n x n covariate
  double g(double y) const { return gw ? gMarGW(y, alpha) : gMarAC(y, alpha); }
  double a(const Chain& ch, int i, int j) const { return dyad[ch.idx(i, j)]; }
 public:
  Reciprocity(bool isGW, Weight w, double al, std::vector<double> d)
    : gw(isGW), weight(w), alpha(al), dyad(std::move(d)) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (i == j) return 0;
    if (weight == NONE) {
      double nIncoming = ch.v(j, i);
      double cont = update * g(nIncoming);
      double nOutgoing = ch.v(i, j);
      return cont + nIncoming * (g(nOutgoing + update) - g(nOutgoing));
    }
    if (weight == DYAD) {
      double nIncoming = ch.v(j, i);
      double cont = update * g(nIncoming) * a(ch, i, j);
      double nOutgoing = ch.v(i, j);
      return cont + nIncoming * (g(nOutgoing + update) * a(ch, j, i) -
                                 g(nOutgoing) * a(ch, j, i));
    }
    double nIncoming = ch.v(j, i) * a(ch, j, i);
    double cont = update * g(nIncoming) * a(ch, i, j);
    double nOutgoing = ch.v(i, j) * a(ch, i, j);
    return cont + nIncoming * (g(nOutgoing + (update * a(ch, i, j))) - g(nOutgoing));
  }
  double target(const Chain& ch, int i, int j) const override {
    if (i == j) return 0;
    if (weight == NONE) return ch.v(i, j) * g(ch.v(j, i));
    if (weight == DYAD) return ch.v(i, j) * g(ch.v(j, i)) * a(ch, i, j);
    double nOutgoing = ch.v(i, j) * a(ch, i, j);
    double nIncoming = ch.v(j, i) * a(ch, j, i);
    return nOutgoing * g(nIncoming);
  }
};

// in_weights_AC, in_weights_GW (R/2_effectFunctions.R): the AC or GW function
// of the inflow of the destination
class InWeights : public Effect {
  bool gw;
  double alpha;
  // sum(cache$valuedNetwork[, j]) as in R (long double accumulation)
  double inWeight(const Chain& ch, int j) const {
    long double total = 0;
    for (int i = 0; i < ch.n; ++i) total += ch.v(i, j);
    return (double)total;
  }
 public:
  InWeights(bool isGW, double a) : gw(isGW), alpha(a) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    double w = inWeight(ch, j);
    double y = update < 0 ? w + update : w;
    return update * (gw ? gMarGW(y, alpha) : gMarAC(y, alpha));
  }
  double target(const Chain& ch, int i, int j) const override {
    double w = inWeight(ch, j);
    return (gw ? gCumGW(w, alpha) : gCumAC(w, alpha)) / ch.n;
  }
};

// ---- transitivity and triad effects (R/2e_effectFunctions_transitivity.R)

// rows i, j and columns i, j of a network, with the loop entries set to 0
struct RowsCols {
  std::vector<double> ri, rj, ci, cj;
  RowsCols(const Chain& ch, const std::vector<double>& m, int i, int j)
    : ri(ch.n), rj(ch.n), ci(ch.n), cj(ch.n) {
    for (int k = 0; k < ch.n; ++k) {
      ri[k] = m[ch.idx(i, k)];
      rj[k] = m[ch.idx(j, k)];
      ci[k] = m[ch.idx(k, i)];
      cj[k] = m[ch.idx(k, j)];
    }
    ri[i] = 0;
    rj[j] = 0;
    ci[i] = 0;
    cj[j] = 0;
  }
};

// R's sum(a * b) and sum(a * b * c): element-wise products, long double sum
inline double sumProd(const std::vector<double>& a, const std::vector<double>& b) {
  long double total = 0;
  for (size_t k = 0; k < a.size(); ++k) total += a[k] * b[k];
  return (double)total;
}
inline double sumProd(const std::vector<double>& a, const std::vector<double>& b,
                      const std::vector<double>& c) {
  long double total = 0;
  for (size_t k = 0; k < a.size(); ++k) total += (a[k] * b[k]) * c[k];
  return (double)total;
}

// transitivity_basic
class TransitivityBasic : public Effect {
 public:
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (i == j) return 0;
    RowsCols n(ch, ch.valued, i, j);
    return update * (sumProd(n.ri, n.cj) + sumProd(n.ri, n.rj) + sumProd(n.ci, n.cj));
  }
  double target(const Chain& ch, int i, int j) const override {
    if (i == j) return 0;
    RowsCols n(ch, ch.valued, i, j);
    return ch.v(i, j) * sumProd(n.ri, n.cj);
  }
};

// transitivity_min
class TransitivityMin : public Effect {
 public:
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (i == j) return 0;
    if (update == 0) return 0;
    double nf = ch.netFlows[ch.idx(i, j)];
    if (update > 0 && nf >= 0) return 0;
    if (update < 0 && nf > 0) return 0;
    RowsCols m(ch, ch.minNet, i, j);
    double dyadValue = ch.minNet[ch.idx(i, j)];
    long double v = 0;
    for (int k = 0; k < ch.n; ++k) {
      double twoPath = std::min(m.ri[k], m.cj[k]);
      if (update > 0) v += dyadValue < twoPath;
      if (update < 0) v += dyadValue <= twoPath;
    }
    return update > 0 ? (double)v : -(double)v;
  }
  double target(const Chain& ch, int i, int j) const override {
    if (i == j) return 0;
    RowsCols m(ch, ch.minNet, i, j);
    double dyadValue = ch.minNet[ch.idx(i, j)];
    long double v = 0;
    for (int k = 0; k < ch.n; ++k) v += std::min(dyadValue, std::min(m.ri[k], m.cj[k]));
    return (double)v / 6;
  }
};

// transitivity_netflow
class TransitivityNetflow : public Effect {
  // sum(pmin(getPathWeights(ik, kj), x)) over k, for row a and column b of
  // the net flow network
  static double sumPathWeights(const Chain& ch, int a, int b, double x) {
    long double total = 0;
    for (int k = 0; k < ch.n; ++k) {
      double ik = ch.netFlows[ch.idx(a, k)], kj = ch.netFlows[ch.idx(k, b)];
      double w = std::min(std::abs(ik), std::abs(kj));
      if (ik < 0 && kj < 0) w = 0;
      total += std::min(w, x);
    }
    return (double)total;
  }
 public:
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (i == j) return 0;
    double nf = ch.netFlows[ch.idx(i, j)];
    double contributionBefore = 0;
    // the path weights are those of i -> k -> j if the net flow is (or
    // becomes) positive, and of j -> k -> i otherwise
    bool forward;
    if (nf > 0) {
      forward = true;
      contributionBefore = sumPathWeights(ch, i, j, nf);
    } else if (nf < 0) {
      forward = false;
      contributionBefore = sumPathWeights(ch, j, i, ch.netFlows[ch.idx(j, i)]);
    } else {
      forward = update > 0;
    }
    double tieWeightAfter = nf + update;
    double contributionAfter = forward ? sumPathWeights(ch, i, j, std::abs(tieWeightAfter))
                                       : sumPathWeights(ch, j, i, std::abs(tieWeightAfter));
    return contributionAfter - contributionBefore;
  }
  double target(const Chain& ch, int i, int j) const override {
    if (i == j) return 0;
    double nf = ch.netFlows[ch.idx(i, j)];
    if (nf <= 0) return 0;
    return sumPathWeights(ch, i, j, nf) / 3;
  }
};

// transitivity_AC, transitivity_GW
class TransitivityACGW : public Effect {
  bool gw;
  double alpha;
  double g(double y) const {
    if (!gw) return 1 - rpow(1 - 1 / alpha, y);
    if (alpha == 1) return y;
    return (1 - rpow(alpha, -y)) / (alpha - 1);
  }
  // R's sum(w * g_mar(y))
  double sumWeighted(const std::vector<double>& w, const std::vector<double>& y) const {
    long double total = 0;
    for (size_t k = 0; k < w.size(); ++k) total += w[k] * g(y[k]);
    return (double)total;
  }
 public:
  TransitivityACGW(bool isGW, double a) : gw(isGW), alpha(a) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (i == j) return 0;
    RowsCols n(ch, ch.valued, i, j);
    int N = ch.n;

    // part 1: i,j closes two-path
    double cont = update * g(sumProd(n.ri, n.cj));

    // part 2: i,j closes instar
    std::vector<double> riU = n.ri;
    riU[j] = riU[j] + update;
    std::vector<double> before(N), after(N);
    for (int k = 0; k < N; ++k) {
      double product = 0;
      for (int m = 0; m < N; ++m) product += n.ri[m] * ch.v(m, k);
      before[k] = product - n.ri[k] * ch.v(k, k);
      after[k] = before[k] + update * n.rj[k];
    }
    double conBe = sumWeighted(n.ri, before) - (n.ri[j] * g(before[j]));
    double conAf = sumWeighted(riU, after) - (riU[j] * g(after[j]));
    cont = cont + (conAf - conBe);

    // part 3: i,j closes outstar
    for (int k = 0; k < N; ++k) {
      double product = 0;
      for (int m = 0; m < N; ++m) product += ch.v(k, m) * n.cj[m];
      before[k] = product - ch.v(k, k) * n.cj[k];
      after[k] = before[k] + update * n.ci[k];
    }
    conBe = sumWeighted(n.cj, before);
    conAf = sumWeighted(n.cj, after);
    cont = cont + (conAf - conBe);
    return cont;
  }
  double target(const Chain& ch, int i, int j) const override {
    if (i == j) return 0;
    RowsCols n(ch, ch.valued, i, j);
    return n.ri[j] * g(sumProd(n.ri, n.cj));
  }
};

// triad120C, triad120D, triad120U
class Triad120 : public Effect {
 public:
  enum Type { C, D, U };
 private:
  Type type;
 public:
  explicit Triad120(Type t) : type(t) {}
  double change(const Chain& ch, int i, int j, int edge, int update) const override {
    if (i == j) return 0;
    RowsCols n(ch, ch.valued, i, j);
    double nji = ch.v(j, i);
    if (type == C) {
      double cont = update * sumProd(n.ri, n.cj) * nji;
      cont = cont + update * sumProd(n.ci, n.rj) * nji;
      cont = cont + update * sumProd(n.rj, n.ri, n.ci);
      return cont + update * sumProd(n.ci, n.cj, n.rj);
    }
    if (type == D) {
      double cont = update * sumProd(n.ci, n.cj) * nji;
      return cont + update * sumProd(n.ri, n.cj, n.rj);
    }
    double cont = update * sumProd(n.ri, n.rj) * nji;
    return cont + update * sumProd(n.cj, n.ri, n.ci);
  }
  double target(const Chain& ch, int i, int j) const override {
    if (i == j) return 0;
    RowsCols n(ch, ch.valued, i, j);
    double recip = ch.v(i, j) * ch.v(j, i);
    if (type == C) return recip * sumProd(n.ri, n.cj);
    if (type == D) return recip * sumProd(n.ci, n.cj) / 2;
    return recip * sumProd(n.ri, n.rj) / 2;
  }
};

}  // namespace

std::unique_ptr<Effect> makeEffect(const Rcpp::List& spec, const Chain& chain) {
  std::string name = Rcpp::as<std::string>(spec["name"]);
  typedef std::unique_ptr<Effect> E;
  if (name == "loops") return E(new Loops());
  if (name == "concentration_AC") {
    return E(new ConcentrationAC(alphaOf(spec, "concentration_AC")));
  }
  if (name == "reciprocity_AC") {
    return E(new ReciprocityAC(alphaOf(spec, "reciprocity_AC")));
  }
  if (name == "dyadic_covariate") return E(new DyadicCovariate(attributeOf(spec)));
  if (name == "alter_covariate") {
    return E(new AlterCovariate(attributeOf(spec), loopContributionOf(spec)));
  }
  if (name == "resource_covar_to_node_covar") {
    return E(new ResourceCovarToNodeCovar(
      chain.resourceIndex(Rcpp::as<std::string>(spec["resourceAttribute"])),
      attributeOf(spec), loopContributionOf(spec)));
  }
  if (name == "loops_resource_covar") {
    return E(new LoopsResourceCovar(
      chain.resourceIndex(Rcpp::as<std::string>(spec["resourceAttribute"]))));
  }
  if (name == "loops_AC" || name == "loops_GW") {
    double alpha = Rcpp::as<double>(spec["alpha"]);
    if (alpha <= 0) Rcpp::stop("Alpha parameter in GW loops weights function must be positive");
    if (name == "loops_AC") return E(new LoopsAC(alpha));
    return E(new LoopsGW(alpha));
  }
  if (name == "loops_additional_origin" || name == "loops_x_loops_additional_origin") {
    return E(new LoopsAdditionalOrigin(
      chain.resourceIndex(Rcpp::as<std::string>(spec["resourceAttribute"])),
      name == "loops_x_loops_additional_origin"));
  }
  if (name == "loops_node_covar") return E(new LoopsNodeCovar(attributeOf(spec)));
  if (name == "loops_resource_covar_node_covar") {
    return E(new LoopsResourceCovarNodeCovar(
      chain.resourceIndex(Rcpp::as<std::string>(spec["resourceAttribute"])),
      attributeOf(spec)));
  }
  if (name == "concentration_GW" || name == "concentration_AC_dyad_covar" ||
      name == "concentration_GW_dyad_covar" ||
      name == "concentration_AC_resource_covar_bin" ||
      name == "concentration_GW_resource_covar_bin") {
    bool gw = name.find("_GW") != std::string::npos;
    Concentration::Weight w = Concentration::NONE;
    std::vector<double> dyad;
    int resource = -1;
    // the R versions name the GW dyadic and unweighted effects
    // "concentration_GW" in their alpha error message
    std::string alphaName = name == "concentration_GW_dyad_covar" ? "concentration_GW" : name;
    double alpha = alphaOf(spec, alphaName.c_str());
    if (name.find("dyad_covar") != std::string::npos) {
      w = Concentration::DYAD;
      dyad = attributeOf(spec);
      if (dyad.size() != (size_t)chain.n * chain.n) {
        Rcpp::stop("attribute.index in %s_bin function must be a square matrix", name);
      }
    }
    if (name.find("resource_covar") != std::string::npos) {
      w = Concentration::RESOURCE;
      resource = chain.resourceIndex(Rcpp::as<std::string>(spec["resourceAttribute"]));
    }
    return E(new Concentration(gw, w, alpha, dyad, resource));
  }
  if (name == "reciprocity_GW" || name == "reciprocity_AC_dyad_covar" ||
      name == "reciprocity_AC_dyad_covar_bin" || name == "reciprocity_GW_dyad_covar" ||
      name == "reciprocity_GW_dyad_covar_bin") {
    bool gw = name.find("_GW") != std::string::npos;
    double alpha = alphaOf(spec, name.c_str());
    Reciprocity::Weight w = Reciprocity::NONE;
    std::vector<double> dyad;
    if (name.find("dyad_covar") != std::string::npos) {
      w = name.find("_bin") != std::string::npos ? Reciprocity::DYAD_BIN : Reciprocity::DYAD;
      dyad = attributeOf(spec);
      if (dyad.size() != (size_t)chain.n * chain.n) {
        Rcpp::stop("attribute.index in %s function must be a square matrix", name);
      }
    }
    return E(new Reciprocity(gw, w, alpha, dyad));
  }
  if (name == "in_weights_AC" || name == "in_weights_GW") {
    double alpha = Rcpp::as<double>(spec["alpha"]);
    if (alpha <= 0) {
      Rcpp::stop("Alpha parameter in %s weights function must be positive", name);
    }
    return E(new InWeights(name == "in_weights_GW", alpha));
  }
  if (name == "transitivity_basic") return E(new TransitivityBasic());
  if (name == "transitivity_min") return E(new TransitivityMin());
  if (name == "transitivity_netflow") return E(new TransitivityNetflow());
  if (name == "transitivity_AC" || name == "transitivity_GW") {
    return E(new TransitivityACGW(name == "transitivity_GW", Rcpp::as<double>(spec["alpha"])));
  }
  if (name == "triad120C") return E(new Triad120(Triad120::C));
  if (name == "triad120D") return E(new Triad120(Triad120::D));
  if (name == "triad120U") return E(new Triad120(Triad120::U));
  Rcpp::stop("effect '%s' has no C++ version", name);
}
