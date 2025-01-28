#pragma once
#ifndef ISLAND_DRIFT
#define ISLAND_DRIFT

#include "component.hxx"

/// Set ion density to a fixed value
///
struct IslandDrift : public Component {
  /// Inputs
  /// - <name>
  ///   - AA
  ///   - Drift velocity
  ///   - density sink proportionality constant (beta in beta*N)
  ///   - momentum diffusion constant (alpha in the term alpha*dN/dx)

  IslandDrift(std::string name, Options& alloptions, Solver *);
  void transform(Options &state) override;
  //void outputVars(Options &state) override;
  
  bool diagnose;
private:
  std::string name; ///< Short name of species e.g "e"

  BoutReal AA;          ///< Atomic mass e.g. proton = 1

  BoutReal vD;         ///< Drift velocity for 1D drift model
  Field3D  sink_propto_N; ///< density sink proportional to N
  BoutReal NV_propto_dN_dx; ///< term proportional to dN/dx in the momentum equation"
};

namespace {
RegisterComponent<IslandDrift> registercomponentIslandDrift("island_drift");
}

#endif // FIXED_DENSITY_H
