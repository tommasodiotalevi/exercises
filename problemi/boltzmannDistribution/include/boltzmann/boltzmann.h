#ifndef BOLTZMANN_H
#define BOLTZMANN_H

#include <random>
#include <string>
#include <vector>

struct SimStats {
  double mean;
  double stddev;
};

class BoltzmannSimulator {
 public:
  BoltzmannSimulator(int molecules, double kBT = 2.0, double hparam = 10.0,
                     unsigned seed = 0);

  void run(int nInteractions);
  const std::vector<double>& getEnergies() const { return energies_; }
  SimStats getStats() const;

 private:
  int molecules_;
  double deltaU_;
  std::vector<double> energies_;

  std::mt19937_64 rng_;
  std::uniform_int_distribution<int> pick_;

  void transferStep();
};

#endif