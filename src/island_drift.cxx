
#include "../include/island_drift.hxx"
#include <bout/fv_ops.hxx>
#include <bout/difops.hxx>
#include "../include/div_ops.hxx"
#include "../include/hermes_build_config.hxx"
#include <bout/constants.hxx>

using bout::globals::mesh;

IslandDrift::IslandDrift(std::string name, Options &alloptions,
                                     Solver *UNUSED(solver)) 
        :name(name){

    auto& options = alloptions[name];

    const Options& units = alloptions["units"];
    const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();

    vD = options["vD"]
      .doc("Drift velocity for 1D models")
      .withDefault(vD);
    sink_propto_N = options["sink_propto_N"]
      .doc("Sink term in ddt(N" + name + std::string(") proportional to N. Units [1/s]"))
      .withDefault(sink_propto_N)
    /Omega_ci;
    NV_propto_dN_dx = options["NV_propto_dN_dx"].doc("term proportional to dN/dx in the momentum equation").withDefault(0.0);

    Theta = options["Theta"]
      .doc("Field-line Pitch defined by Feng et al.")
      .withDefault(1e-3);
    Theta_inv = 1./Theta;
    
    diagnose = options["diagnose"]
      .doc("Save additional output diagnostics")
      .withDefault<bool>(false);
}

void IslandDrift::transform(Options &state) {
    AUTO_TRACE();
    
    auto& species = state["species"][name];
    auto AA = get<BoutReal>(species["AA"]);

    const Field3D N = get<Field3D>(species["density"]);

    // Typical wave speed used for numerical diffusion
    Field3D fastest_wave;
    if (state.isSet("fastest_wave")) {
        fastest_wave = get<Field3D>(state["fastest_wave"]);
    } else {
        Field3D T = get<Field3D>(species["temperature"]);
        fastest_wave = sqrt(T / AA);
    }

    Field3D density_drift(0.0);
    Field3D momentum_drift(0.0);
    //Field3D energy_drift(0.0);

    density_drift  -= Theta_inv*Grad_par(vD*N);

    if (species.isSet("velocity")) {
        Field3D V = get<Field3D>(species["velocity"]);
        Field3D NV = get<Field3D>(species["momentum"]);
        momentum_drift -= AA*Theta_inv*FV::Div_par_mod<hermes::Limiter>(vD*N, V, fastest_wave, flow_ylow);
    } //add a check for 1D later

    density_drift -= sink_propto_N * N;
    momentum_drift -= AA*NV_propto_dN_dx*Grad_par(N);

    add(species["density_source"], density_drift);
    add(species["momentum_source"], momentum_drift); //should check if evolve momentum is on.
    //add(species["energy_source"], energy_drift);
}

/*
void IslandDrift::outputVars(Options& state) {
  AUTO_TRACE();

  if (diagnose) {
    auto Nnorm = get<BoutReal>(state["Nnorm"]);
    auto Tnorm = get<BoutReal>(state["Tnorm"]);
    auto Omega_ci = get<BoutReal>(state["Omega_ci"]);

    BoutReal Cs0 = sqrt(SI::qe * Tnorm / (AA * Mp));
    BoutReal DivQnorm = SI::qe * Tnorm * Nnorm * Omega_ci;

    set_with_attrs(state["vD_grad_par_N"], Theta_inv*vD*Grad_par(N),
                   {{"time_dimension", "t"},
                    {"units", "1 / (m^3 s)"},
                    {"conversion", Nnorm*Omega_ci},
                    {"long_name", "Drift velocity"},
                    {"source", "island_drift"}});

    set_with_attrs(state["vD_Div_par_NV"], sink_propto_N * N;,
                   {{"time_dimension", "t"},
                    {"units", "1 / (m^3 s)},
                    {"conversion", Nnorm*Omega_ci},
                    {"long_name", "Private flux region sink"},
                    {"source", "island_drift"}});

    set_with_attrs(state["vD_Div_par_NV"], NV_propto_dN_dx*Grad_par(N),
                   {{"time_dimension", "t"},
                    {"units", "1 / (m^2 s^2)"},
                    {"conversion", Cs0*Nnorm*Omega_ci},
                    {"long_name", "Momentum term proportional to dn/dx"},
                    {"source", "island_drift"}});


  }
}
*/