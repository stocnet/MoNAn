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

}  // namespace

std::unique_ptr<Effect> makeEffect(const Rcpp::List& spec, const Chain& chain) {
  std::string name = Rcpp::as<std::string>(spec["name"]);
  if (name == "loops") return std::unique_ptr<Effect>(new Loops());
  Rcpp::stop("effect '%s' has no C++ version", name);
}
