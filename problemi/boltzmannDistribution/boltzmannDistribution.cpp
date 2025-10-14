#include <TCanvas.h>
#include <TF1.h>
#include <TH1D.h>
#include <TLegend.h>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "boltzmann/boltzmann.h"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

void plotDistributionROOT(const std::vector<int> &interactions_list,
                          int molecules) {
  TCanvas *c1 = new TCanvas("c1", "Boltzmann Simulations", 1200, 800);
  c1->Divide(std::ceil(interactions_list.size() / 2.), 2);
  std::cout << "aa" << std::ceil(interactions_list.size() / 2) << std::endl;
  int pad{1};
  for (auto nInteractions : interactions_list) {
    BoltzmannSimulator sim(molecules);
    sim.run(nInteractions);

    SimStats st{sim.getStats()};
    std::cout << "Number of interactions=" << nInteractions
              << " mean=" << st.mean << " stddev=" << st.stddev << std::endl;

    c1->cd(pad);
    std::string hname{"hE_" + std::to_string(nInteractions)};
    TH1D *hE = new TH1D(
        hname.c_str(),
        ("Energy distribution interactions=" + std::to_string(nInteractions))
            .c_str(),
        80, 0, 10);
    for (double e : sim.getEnergies()) hE->Fill(e);

    hE->GetXaxis()->SetTitle("Energy");
    hE->GetYaxis()->SetTitle("Count");
    hE->SetLineColor(pad);
    hE->Draw();

    // Perform the exponential fit
    if (nInteractions > molecules * 10) {
      TF1 *fexp = new TF1("fexp", "expo", 0, 10);
      hE->Fit(fexp, "RQ");
      fexp->SetLineColor(2);
      fexp->Draw("same");
    }
    pad++;
  }

  c1->SaveAs("problemi/boltzmannDistribution/boltzmann_distributions.png");
}

void plotDistributionSTDOUT(const std::vector<int> &interaction_list,
                            int molecules) {
  for (auto nInteractions : interaction_list) {
    BoltzmannSimulator sim(molecules);
    sim.run(nInteractions);

    SimStats st{sim.getStats()};
    std::cout << "\n=== Interaction count: " << nInteractions
              << " ===" << std::endl;
    std::cout << "Mean energy:   " << st.mean << std::endl;
    std::cout << "Std deviation: " << st.stddev << std::endl;

    const int nbins{20};
    const double emin{0.0}, emax{10.0};
    std::vector<int> counts(nbins, 0);

    for (double e : sim.getEnergies()) {
      if (e >= emin && e < emax) {
        int bin = static_cast<int>((e - emin) / (emax - emin) * nbins);
        if (bin >= 0 && bin < nbins) counts[bin]++;
      }
    }

    std::cout << "Energy distribution:\n";
    for (int i = 0; i < nbins; ++i) {
      double low{emin + i * (emax - emin) / nbins};
      double high{low + (emax - emin) / nbins};
      int barlen = static_cast<int>(
          50.0 * counts[i] /
          (*std::max_element(counts.begin(), counts.end()) + 1e-9));

      std::cout << std::fixed << std::setprecision(1) << "[" << std::setw(4)
                << low << " - " << std::setw(4) << high << "] "
                << std::string(barlen, '#') << " (" << counts[i] << ")\n";
    }
  }
}

int boltzmannDistribution(int argc, char **argv) {
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0]
              << " molecules interactions[list(comma-separated)]"
              << "Example: " << argv[0] << " 1000 1000,10000,100000"
              << std::endl;
    return 1;
  }

  int molecules{std::stoi(argv[1])};
  std::string s_interaction_list{argv[2]};

  std::vector<int> interaction_list;
  std::stringstream ss(s_interaction_list);
  std::string token;
  while (std::getline(ss, token, ',')) {
    if (!token.empty())
      interaction_list.push_back(
          std::stoll(token));  // Parse from list of interactions (e.g.
                               // "1000,10000,100000")
  }

  // Plot using ROOT
  // plotDistributionROOT(interaction_list, molecules);

  // Plot using STDOUT
  plotDistributionSTDOUT(interaction_list, molecules);

  return 0;
}

TEST_CASE(
    "Test total energy conservation and that energies are always positive") {
  int molecules{1000};
  int nInteractions{100000};
  double kBT{2.0};

  double initialEnergy{molecules * 1.5 * kBT};
  BoltzmannSimulator sim(molecules, kBT);
  sim.run(nInteractions);

  double finalEnergy{0.0};

  for (double e : sim.getEnergies()) {
    finalEnergy += e;
    CHECK(e >= 0.0);
  }
  CHECK(finalEnergy == doctest::Approx(initialEnergy).epsilon(1e-9));
}