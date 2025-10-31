#include "boltzmann/boltzmann.h"

#include <algorithm>
#include <cmath>

BoltzmannSimulator::BoltzmannSimulator(int molecules, double kBT, double hparam,
                                       unsigned seed)
    : molecules_{molecules},
      deltaU_{kBT / hparam},
      energies_(molecules, 1.5 * kBT),
      rng_{seed ? seed : std::random_device{}()},
      pick_{0, molecules - 1} {}

void BoltzmannSimulator::transferStep() {
  int i{pick_(rng_)};
  int j{pick_(rng_)};
  if (i == j) return;  // For the moment, we exclude self-transfers
  if (energies_[i] >= deltaU_) {
    energies_[i] -= deltaU_;
    energies_[j] += deltaU_;
  }
}

void BoltzmannSimulator::run(int nInteractions) {
  for (int it = 0; it < nInteractions; ++it) {
    transferStep();
  }
}

SimStats BoltzmannSimulator::getStats() const {
  SimStats s{0., 0.};
  double sum{0.};
  for (double e : energies_) sum += e;
  s.mean = sum / molecules_;
  double var{0.};
  for (double e : energies_) var += (e - s.mean) * (e - s.mean);
  s.stddev = std::sqrt(var / molecules_);
  return s;
}