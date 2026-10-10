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
  Rcpp::stop("effect '%s' has no C++ version", name);
}
