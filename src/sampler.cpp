// C++ implementation of the MCMC sampler (see R/1b_cppSampler.R).
// One step follows simulateOneStep() in R/1_hiddenFunctions.R exactly,
// including the order of random numbers, so that for the same seed the C++
// chain equals the R chain.
#include "monan.h"

#include <algorithm>
#include <map>

namespace {

// R's sample(x, 1): a single number x >= 1 means sample(1:x, 1)
int sampleLikeR(const std::vector<int>& x) {
  if (x.size() == 1 && x[0] >= 1) return (int)R_unif_index((double)x[0]) + 1;
  return x[(size_t)R_unif_index((double)x.size())];
}

std::vector<double> asVector(SEXP m) {
  Rcpp::NumericVector v(m);
  return std::vector<double>(v.begin(), v.end());
}

Rcpp::NumericMatrix asMatrix(const std::vector<double>& v, int n) {
  Rcpp::NumericMatrix m(n, n);
  std::copy(v.begin(), v.end(), m.begin());
  return m;
}

class Sampler {
 public:
  Chain ch;
  std::vector<std::unique_ptr<Effect>> effects;
  std::vector<double> parameters;
  std::vector<int> resourceSet, nodeSet;
  bool allowLoops, multinomial;

  void step() {
    int edge = sampleLikeR(resourceSet) - 1;
    int i = ch.sender[edge], j = ch.receiver[edge];

    std::vector<int> nodes;
    if (allowLoops) {
      nodes = nodeSet;
    } else {
      for (int k : nodeSet) if (k != i + 1) nodes.push_back(k);
    }

    size_t K = effects.size();

    if (multinomial) {
      multinomialStep(edge, i, j, nodes);
      return;
    }

    int l = sampleLikeR(nodes) - 1;
    if (j == l) return;

    std::vector<double> drop(K), create(K);
    for (size_t k = 0; k < K; ++k) drop[k] = effects[k]->change(ch, i, j, edge, -1);

    // remove the resource, keeping what is needed to undo this if rejected
    Chain::Cells old = ch.getCells(i, j);
    ch.update(i, j, edge, -1);
    ch.sender[edge] = ch.receiver[edge] = -1;

    for (size_t k = 0; k < K; ++k) create[k] = effects[k]->change(ch, i, l, edge, 1);

    long double contribution = 0;
    for (size_t k = 0; k < K; ++k) contribution += (drop[k] + create[k]) * parameters[k];
    double acceptance = std::min(1.0, std::exp((double)contribution));

    if (acceptance >= unif_rand()) {
      ch.update(i, l, edge, 1);
      ch.sender[edge] = i;
      ch.receiver[edge] = l;
    } else {
      ch.restoreCells(i, j, old);
      ch.sender[edge] = i;
      ch.receiver[edge] = j;
    }
  }

  // multinomial proposal: the destination is drawn from all possible ones with
  // probabilities proportional to exp(parameters . change statistics), with the
  // same arithmetic as R (colSums, sum and cumsum accumulate in long double)
  void multinomialStep(int edge, int i, int j, const std::vector<int>& nodes) {
    size_t K = effects.size();
    std::vector<double> drop(K);
    for (size_t k = 0; k < K; ++k) drop[k] = effects[k]->change(ch, i, j, edge, -1);

    ch.update(i, j, edge, -1);
    ch.sender[edge] = ch.receiver[edge] = -1;

    size_t D = nodes.size();
    std::vector<double> expValues(D);
    long double total = 0;
    for (size_t d = 0; d < D; ++d) {
      long double objective = 0;
      for (size_t k = 0; k < K; ++k) {
        double create = effects[k]->change(ch, i, nodes[d] - 1, edge, 1);
        objective += (create + drop[k]) * parameters[k];
      }
      expValues[d] = std::exp((double)objective);
      total += expValues[d];
    }

    double u = unif_rand();
    long double cumulative = 0;
    size_t draw = D;
    for (size_t d = 0; d < D; ++d) {
      cumulative += expValues[d] / (double)total;
      if ((double)cumulative >= u) {
        draw = d;
        break;
      }
    }
    if (draw == D) Rcpp::stop("no destination drawn in multinomial proposal");
    int l = nodes[draw] - 1;

    ch.update(i, l, edge, 1);
    ch.sender[edge] = i;
    ch.receiver[edge] = l;
  }

  // as getNetworkStatistics(): sum of target contributions over all cells,
  // in the order of expand.grid (i fastest)
  std::vector<double> statistics() const {
    std::vector<double> s(effects.size());
    for (size_t k = 0; k < effects.size(); ++k) s[k] = effects[k]->targetTotal(ch);
    return s;
  }

  Rcpp::IntegerMatrix edgelist() const {
    Rcpp::IntegerMatrix m(ch.sender.size(), 2);
    for (size_t e = 0; e < ch.sender.size(); ++e) {
      m(e, 0) = ch.sender[e] + 1;
      m(e, 1) = ch.receiver[e] + 1;
    }
    return m;
  }

  Rcpp::List cache() const {
    Rcpp::List resource(ch.resourceNames.size());
    for (size_t k = 0; k < ch.resourceNames.size(); ++k) {
      resource[k] = asMatrix(ch.resourceNetworks[k], ch.n);
    }
    resource.attr("names") = ch.resourceNames;
    return Rcpp::List::create(
      Rcpp::Named("valuedNetwork") = asMatrix(ch.valued, ch.n),
      Rcpp::Named("resourceNetworks") = resource,
      Rcpp::Named("netFlowsNetwork") = asMatrix(ch.netFlows, ch.n),
      Rcpp::Named("minNetwork") = asMatrix(ch.minNet, ch.n),
      Rcpp::Named("twoPathNetwork") = asMatrix(ch.twoPaths, ch.n));
  }
};

// sets up chain and effects from the R objects
void setupSampler(Sampler& s,
                  const Rcpp::IntegerMatrix& edgelist,
                  const Rcpp::List& cache,
                  const Rcpp::List& resourceCovariates,
                  const Rcpp::List& effectSpecs) {
  Rcpp::NumericMatrix valued = cache["valuedNetwork"];
  s.ch.n = valued.nrow();
  for (int e = 0; e < edgelist.nrow(); ++e) {
    s.ch.sender.push_back(edgelist(e, 0) - 1);
    s.ch.receiver.push_back(edgelist(e, 1) - 1);
  }
  s.ch.valued = asVector(cache["valuedNetwork"]);
  s.ch.netFlows = asVector(cache["netFlowsNetwork"]);
  s.ch.minNet = asVector(cache["minNetwork"]);
  s.ch.twoPaths = asVector(cache["twoPathNetwork"]);
  if (cache.containsElementNamed("resourceNetworks")) {
    Rcpp::List rn = cache["resourceNetworks"];
    if (rn.size() > 0) {
      Rcpp::CharacterVector names = rn.names();
      for (int k = 0; k < rn.size(); ++k) {
        std::string name = Rcpp::as<std::string>(names[k]);
        s.ch.resourceNames.push_back(name);
        s.ch.resourceNetworks.push_back(asVector(rn[k]));
        s.ch.resourceValues.push_back(asVector(resourceCovariates[name]));
      }
    }
  }
  for (int k = 0; k < effectSpecs.size(); ++k) {
    s.effects.push_back(makeEffect(Rcpp::as<Rcpp::List>(effectSpecs[k]), s.ch));
  }
}

}  // namespace

// Runs one chain: burnIn steps, then `iterations` times (thinning + 1) steps
// followed by the calculation of the statistics, as simulateStatisticVectors().
// With iterations = 0, this is simulateNSteps(n = burnIn).
// [[Rcpp::export]]
Rcpp::List simulateChainCpp(Rcpp::IntegerMatrix edgelist,
                            Rcpp::List cache,
                            Rcpp::List resourceCovariates,
                            Rcpp::IntegerVector resourceSet,
                            Rcpp::IntegerVector nodeSet,
                            Rcpp::List effectSpecs,
                            Rcpp::NumericVector parameters,
                            int burnIn,
                            int iterations,
                            int thinning,
                            bool allowLoops,
                            bool multinomialProposal,
                            bool returnDraws) {
  Sampler s;
  setupSampler(s, edgelist, cache, resourceCovariates, effectSpecs);
  s.parameters.assign(parameters.begin(), parameters.end());
  if (s.parameters.size() != s.effects.size()) {
    Rcpp::stop("number of parameters and effects differ");
  }
  s.resourceSet.assign(resourceSet.begin(), resourceSet.end());
  s.nodeSet.assign(nodeSet.begin(), nodeSet.end());
  s.allowLoops = allowLoops;
  s.multinomial = multinomialProposal;

  for (int t = 0; t < burnIn; ++t) s.step();

  Rcpp::NumericMatrix stats(iterations, (int)s.effects.size());
  Rcpp::List draws(returnDraws ? iterations : 0);
  for (int it = 0; it < iterations; ++it) {
    for (int t = 0; t < thinning + 1; ++t) s.step();
    std::vector<double> st = s.statistics();
    for (size_t k = 0; k < st.size(); ++k) stats(it, k) = st[k];
    if (returnDraws) {
      draws[it] = Rcpp::List::create(Rcpp::Named("edgelist") = s.edgelist(),
                                     Rcpp::Named("cache") = s.cache());
    }
  }

  return Rcpp::List::create(Rcpp::Named("statistics") = stats,
                            Rcpp::Named("edgelist") = s.edgelist(),
                            Rcpp::Named("cache") = s.cache(),
                            Rcpp::Named("draws") = draws);
}

// target statistics of a state, as getNetworkStatistics()
// [[Rcpp::export]]
Rcpp::NumericVector networkStatisticsCpp(Rcpp::IntegerMatrix edgelist,
                                         Rcpp::List cache,
                                         Rcpp::List resourceCovariates,
                                         Rcpp::List effectSpecs) {
  Sampler s;
  setupSampler(s, edgelist, cache, resourceCovariates, effectSpecs);
  std::vector<double> st = s.statistics();
  return Rcpp::NumericVector(st.begin(), st.end());
}

// change statistic and target contribution of one effect, for tests
// [[Rcpp::export]]
Rcpp::NumericVector effectValuesCpp(Rcpp::IntegerMatrix edgelist,
                                    Rcpp::List cache,
                                    Rcpp::List resourceCovariates,
                                    Rcpp::List effectSpec,
                                    int i, int j, int edge, int update) {
  Sampler s;
  setupSampler(s, edgelist, cache, resourceCovariates, Rcpp::List::create(effectSpec));
  return Rcpp::NumericVector::create(
    Rcpp::Named("change") = s.effects[0]->change(s.ch, i - 1, j - 1, edge - 1, update),
    Rcpp::Named("target") = s.effects[0]->target(s.ch, i - 1, j - 1));
}

// sums of a resource covariate per tie, adding the values in the order of
// the edge list starting from 0 (as createInternalCache did in R)
// [[Rcpp::export]]
Rcpp::NumericVector tieSumsCpp(Rcpp::IntegerVector tie, Rcpp::NumericVector values, int nCells) {
  Rcpp::NumericVector sums(nCells);
  for (R_xlen_t e = 0; e < tie.size(); ++e) sums[tie[e] - 1] += values[e];
  return sums;
}

// Log pseudo-likelihood of the conditional logit model of each resource's
// destination given all other resources (as getMPLE() in R), and optionally
// its gradient and Hessian, for parameters theta. For each resource, the
// change statistics of adding it to every possible destination are computed
// after removing it (as in the multinomial proposal); nothing is stored.
// [[Rcpp::export]]
Rcpp::List mpleNewtonCpp(Rcpp::IntegerMatrix edgelist,
                         Rcpp::List cache,
                         Rcpp::List resourceCovariates,
                         Rcpp::IntegerVector resourceSet,
                         Rcpp::IntegerVector nodeSet,
                         Rcpp::List effectSpecs,
                         Rcpp::NumericVector theta,
                         bool allowLoops,
                         bool derivatives) {
  Sampler s;
  setupSampler(s, edgelist, cache, resourceCovariates, effectSpecs);
  size_t K = s.effects.size();
  if ((size_t)theta.size() != K) Rcpp::stop("number of parameters and effects differ");

  long double logPL = 0;
  std::vector<double> gradient(K, 0.0), hessian(K * K, 0.0);
  std::vector<double> X, eta, prob, expected(K);

  for (int r : resourceSet) {
    int edge = r - 1;
    int i = s.ch.sender[edge], j = s.ch.receiver[edge];
    std::vector<int> destinations;
    for (int k : nodeSet) if (allowLoops || k != i + 1) destinations.push_back(k - 1);
    size_t D = destinations.size();
    int observed = -1;
    for (size_t d = 0; d < D; ++d) if (destinations[d] == j) observed = (int)d;
    if (observed < 0) {
      // the observed destination is not possible (as match() giving NA in R)
      return Rcpp::List::create(Rcpp::Named("logPL") = NA_REAL);
    }

    // change statistics of adding the resource to each destination
    Chain::Cells old = s.ch.getCells(i, j);
    s.ch.update(i, j, edge, -1);
    s.ch.sender[edge] = s.ch.receiver[edge] = -1;
    X.assign(D * K, 0.0);
    for (size_t d = 0; d < D; ++d) {
      for (size_t k = 0; k < K; ++k) {
        X[d * K + k] = s.effects[k]->change(s.ch, i, destinations[d], edge, 1);
      }
    }
    s.ch.restoreCells(i, j, old);
    s.ch.sender[edge] = i;
    s.ch.receiver[edge] = j;

    // conditional logit contributions (softmax with the maximum subtracted)
    eta.assign(D, 0.0);
    double maxEta = R_NegInf;
    for (size_t d = 0; d < D; ++d) {
      for (size_t k = 0; k < K; ++k) eta[d] += X[d * K + k] * theta[k];
      maxEta = std::max(maxEta, eta[d]);
    }
    double total = 0;
    prob.assign(D, 0.0);
    for (size_t d = 0; d < D; ++d) {
      prob[d] = std::exp(eta[d] - maxEta);
      total += prob[d];
    }
    logPL += eta[observed] - maxEta - std::log(total);
    if (!derivatives) continue;

    for (size_t d = 0; d < D; ++d) prob[d] /= total;
    std::fill(expected.begin(), expected.end(), 0.0);
    for (size_t d = 0; d < D; ++d) {
      for (size_t k = 0; k < K; ++k) expected[k] += prob[d] * X[d * K + k];
    }
    for (size_t k = 0; k < K; ++k) gradient[k] += X[observed * K + k] - expected[k];
    for (size_t a = 0; a < K; ++a) {
      for (size_t b = 0; b < K; ++b) {
        double second = 0;
        for (size_t d = 0; d < D; ++d) second += prob[d] * X[d * K + a] * X[d * K + b];
        hessian[a * K + b] -= second - expected[a] * expected[b];
      }
    }
  }

  Rcpp::NumericMatrix H(K, K);
  for (size_t a = 0; a < K; ++a) {
    for (size_t b = 0; b < K; ++b) H(a, b) = hessian[a * K + b];
  }
  return Rcpp::List::create(
    Rcpp::Named("logPL") = (double)logPL,
    Rcpp::Named("gradient") = Rcpp::NumericVector(gradient.begin(), gradient.end()),
    Rcpp::Named("hessian") = H);
}

// Choice sets of the conditional logit model of the MPLE (as mpleNewtonCpp),
// computed once and stored, each distinct choice set only once. Resources with
// the same origin, the same observed destination and the same values of the
// resource covariates in `keyCovariates` (those used by the model's effects)
// have the same change statistics, so one of them represents the group, with
// the group's size as weight. Returns possible = FALSE if an observed
// destination is not a possible destination (as match() giving NA in R), and
// tooLarge = TRUE, without computing the change statistics, if they would need
// more than maxBytes.
// [[Rcpp::export]]
Rcpp::List mpleChoiceSetsCpp(Rcpp::IntegerMatrix edgelist,
                             Rcpp::List cache,
                             Rcpp::List resourceCovariates,
                             Rcpp::IntegerVector resourceSet,
                             Rcpp::IntegerVector nodeSet,
                             Rcpp::List effectSpecs,
                             Rcpp::CharacterVector keyCovariates,
                             bool allowLoops,
                             double maxBytes) {
  Sampler s;
  setupSampler(s, edgelist, cache, resourceCovariates, effectSpecs);
  size_t K = s.effects.size();
  std::vector<int> keyIndex;
  for (R_xlen_t k = 0; k < keyCovariates.size(); ++k) {
    keyIndex.push_back(s.ch.resourceIndex(Rcpp::as<std::string>(keyCovariates[k])));
  }

  // group the resources, in the order of resourceSet
  std::map<std::vector<double>, size_t> groupOf;
  std::vector<int> representative;
  std::vector<double> weight;
  for (int r : resourceSet) {
    int edge = r - 1;
    int i = s.ch.sender[edge], j = s.ch.receiver[edge];
    bool possible = (allowLoops || i != j) &&
      std::find(nodeSet.begin(), nodeSet.end(), j + 1) != nodeSet.end();
    if (!possible) {
      return Rcpp::List::create(Rcpp::Named("possible") = false);
    }
    std::vector<double> key = {(double)i, (double)j};
    bool missing = false;
    for (int k : keyIndex) {
      double value = s.ch.resourceValues[k][edge];
      if (ISNAN(value)) missing = true;
      key.push_back(value);
    }
    if (!missing) {
      auto found = groupOf.find(key);
      if (found != groupOf.end()) {
        weight[found->second] += 1;
        continue;
      }
      groupOf[key] = representative.size();
    }
    // a new group (resources with missing values are not grouped)
    representative.push_back(edge);
    weight.push_back(1);
  }

  size_t G = representative.size();
  std::vector<size_t> nDestinations(G);
  double rows = 0;
  for (size_t g = 0; g < G; ++g) {
    int i = s.ch.sender[representative[g]];
    size_t D = 0;
    for (int k : nodeSet) if (allowLoops || k != i + 1) ++D;
    nDestinations[g] = D;
    rows += D;
  }
  if (rows * K * sizeof(double) > maxBytes) {
    return Rcpp::List::create(Rcpp::Named("possible") = true,
                              Rcpp::Named("tooLarge") = true);
  }

  // change statistics of adding each representative to each destination,
  // after removing it; one column per destination (effects in rows)
  Rcpp::NumericMatrix X((int)K, (int)rows);
  Rcpp::IntegerVector start((int)G + 1), observed((int)G);
  int row = 0;
  for (size_t g = 0; g < G; ++g) {
    int edge = representative[g];
    int i = s.ch.sender[edge], j = s.ch.receiver[edge];
    start[g] = row;
    Chain::Cells old = s.ch.getCells(i, j);
    s.ch.update(i, j, edge, -1);
    s.ch.sender[edge] = s.ch.receiver[edge] = -1;
    for (int k : nodeSet) {
      int destination = k - 1;
      if (!allowLoops && destination == i) continue;
      if (destination == j) observed[g] = row;
      for (size_t e = 0; e < K; ++e) {
        X(e, row) = s.effects[e]->change(s.ch, i, destination, edge, 1);
      }
      ++row;
    }
    s.ch.restoreCells(i, j, old);
    s.ch.sender[edge] = i;
    s.ch.receiver[edge] = j;
  }
  start[G] = row;

  return Rcpp::List::create(
    Rcpp::Named("possible") = true,
    Rcpp::Named("tooLarge") = false,
    Rcpp::Named("X") = X,
    Rcpp::Named("start") = start,
    Rcpp::Named("observed") = observed,
    Rcpp::Named("weight") = Rcpp::NumericVector(weight.begin(), weight.end()));
}

// Log pseudo-likelihood, and optionally gradient and Hessian, from the stored
// choice sets of mpleChoiceSetsCpp: as mpleNewtonCpp, with each choice set's
// contribution multiplied by its weight.
// [[Rcpp::export]]
Rcpp::List mpleEvaluateCpp(Rcpp::NumericMatrix X,
                           Rcpp::IntegerVector start,
                           Rcpp::IntegerVector observed,
                           Rcpp::NumericVector weight,
                           Rcpp::NumericVector theta,
                           bool derivatives) {
  size_t K = X.nrow();
  if ((size_t)theta.size() != K) Rcpp::stop("number of parameters and effects differ");
  const double* x = X.begin();
  long double logPL = 0;
  std::vector<double> gradient(K, 0.0), hessian(K * K, 0.0);
  std::vector<double> eta, prob, expected(K);

  for (R_xlen_t g = 0; g < weight.size(); ++g) {
    int first = start[g], D = start[g + 1] - start[g];
    int obs = observed[g] - first;
    double w = weight[g];
    const double* xg = x + (size_t)first * K;

    // conditional logit contributions (softmax with the maximum subtracted)
    eta.assign(D, 0.0);
    double maxEta = R_NegInf;
    for (int d = 0; d < D; ++d) {
      for (size_t k = 0; k < K; ++k) eta[d] += xg[d * K + k] * theta[k];
      maxEta = std::max(maxEta, eta[d]);
    }
    double total = 0;
    prob.assign(D, 0.0);
    for (int d = 0; d < D; ++d) {
      prob[d] = std::exp(eta[d] - maxEta);
      total += prob[d];
    }
    logPL += w * (eta[obs] - maxEta - std::log(total));
    if (!derivatives) continue;

    for (int d = 0; d < D; ++d) prob[d] /= total;
    std::fill(expected.begin(), expected.end(), 0.0);
    for (int d = 0; d < D; ++d) {
      for (size_t k = 0; k < K; ++k) expected[k] += prob[d] * xg[d * K + k];
    }
    for (size_t k = 0; k < K; ++k) gradient[k] += w * (xg[obs * K + k] - expected[k]);
    for (size_t a = 0; a < K; ++a) {
      for (size_t b = 0; b < K; ++b) {
        double second = 0;
        for (int d = 0; d < D; ++d) second += prob[d] * xg[d * K + a] * xg[d * K + b];
        hessian[a * K + b] -= w * (second - expected[a] * expected[b]);
      }
    }
  }

  if (!derivatives) {
    return Rcpp::List::create(Rcpp::Named("logPL") = (double)logPL);
  }
  Rcpp::NumericMatrix H(K, K);
  for (size_t a = 0; a < K; ++a) {
    for (size_t b = 0; b < K; ++b) H(a, b) = hessian[a * K + b];
  }
  return Rcpp::List::create(
    Rcpp::Named("logPL") = (double)logPL,
    Rcpp::Named("gradient") = Rcpp::NumericVector(gradient.begin(), gradient.end()),
    Rcpp::Named("hessian") = H);
}
