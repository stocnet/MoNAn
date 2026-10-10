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
  Rcpp::stop("effect '%s' has no C++ version", name);
}
