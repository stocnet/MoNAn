// C++ implementation of the MCMC sampler (see R/1b_cppSampler.R).
// One step follows simulateOneStep() in R/1_hiddenFunctions.R exactly,
// including the order of random numbers, so that for the same seed the C++
// chain equals the R chain.
#include "monan.h"

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
      Rcpp::Named("minNetwork") = asMatrix(ch.minNet, ch.n));
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
