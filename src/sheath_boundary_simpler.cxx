#include "../include/sheath_boundary_simpler.hxx"

#include "bout/constants.hxx"
#include "bout/mesh.hxx"
using bout::globals::mesh;

namespace {
BoutReal clip(BoutReal value, BoutReal min, BoutReal max) {
  if (value < min)
    return min;
  if (value > max)
    return max;
  return value;
}

BoutReal floor(BoutReal value, BoutReal min) {
  if (value < min)
    return min;
  return value;
}

Ind3D indexAt(const Field3D& f, int x, int y, int z) {
  int ny = f.getNy();
  int nz = f.getNz();
  return Ind3D{(x * ny + y) * nz + z, ny, nz};
}

/// Limited free gradient of log of a quantity
/// This ensures that the guard cell values remain positive
/// while also ensuring that the quantity never increases
///
///  fm  fc | fp
///         ^ boundary
/// 
/// exp( 2*log(fc) - log(fm) )
/// Mode 0: default (exponential extrapolation if decreases, Neumann if increases)
/// Mode 1: always exponential extrapolation
/// Mode 2: always linear extrapolation

BoutReal limitFree(BoutReal fm, BoutReal fc, BoutReal mode) {
  if ((fm < fc) && (mode == 0)) {
    return fc; // Neumann rather than increasing into boundary
  }
  if (fm < 1e-10) {
    return fc; // Low / no density condition
  }

  BoutReal fp = 0;
  if ((mode == 0) || (mode == 1)) {
    fp = SQ(fc) / fm;     // Exponential
  } else if (mode == 2) {
    fp = 2.0 * fc - fm;   // Linear
  } else {
    throw BoutException("Unknown boundary mode");
  }

  return fp;  // Extrapolation


#if CHECKLEVEL >= 2
  if (!std::isfinite(fp)) {
    throw BoutException("SheathBoundary limitFree: {}, {} -> {}", fm, fc, fp);
  }
#endif

  return fp;
}

} // namespace

SheathBoundarySimpler::SheathBoundarySimpler(std::string name, Options& alloptions,
                                           Solver*) {
  AUTO_TRACE();

  Options& options = alloptions[name];

  vD = options["vD"].doc("Drift velocity for 1D models").withDefault(-1.0);

  sheath_ion_polytropic = options["sheath_ion_polytropic"]
         .doc("Ion polytropic coefficient in Bohm sound speed")
         .withDefault(1.0);

  lower_y = options["lower_y"].doc("Boundary on lower y?").withDefault<bool>(true);
  upper_y = options["upper_y"].doc("Boundary on upper y?").withDefault<bool>(true);


  const Options& units = alloptions["units"];
  const BoutReal Tnorm = units["eV"];

  no_flow = options["no_flow"]
    .doc("Set zero particle flow, keeping energy flow")
    .withDefault<bool>(false);

  density_boundary_mode = options["density_boundary_mode"]
    .doc("BC mode: 0=LimitFree, 1=ExponentialFree, 2=LinearFree")
    .withDefault<BoutReal>(1);

  pressure_boundary_mode = options["pressure_boundary_mode"]
    .doc("BC mode: 0=LimitFree, 1=ExponentialFree, 2=LinearFree")
    .withDefault<BoutReal>(1);

  temperature_boundary_mode = options["temperature_boundary_mode"]
    .doc("BC mode: 0=LimitFree, 1=ExponentialFree, 2=LinearFree")
    .withDefault<BoutReal>(1);

  
}

void SheathBoundarySimpler::transform(Options& state) {
  AUTO_TRACE();

  Options& allspecies = state["species"];
  Coordinates* coord = mesh->getCoordinates();

  //////////////////////////////////////////////////////////////////
  // Electrostatic potential
  // If phi is set, use free boundary condition
  // If phi not set, calculate assuming zero current
  Field3D phi;
  if (IS_SET_NOBOUNDARY(state["fields"]["phi"])) {
    phi = toFieldAligned(getNoBoundary<Field3D>(state["fields"]["phi"]));
  } else {
    // Calculate potential phi assuming zero current

    // Need to sum  n_i Z_i C_i over all ion species
    //
    // To avoid looking up species for every grid point, this
    // loops over the boundaries once per species.
    Field3D ion_sum = 0.0;

    // Iterate through charged ion species
    for (auto& kv : allspecies.getChildren()) {
      const Options& species = kv.second;

      if ((kv.first == "e") or !species.isSet("charge")
          or (get<BoutReal>(species["charge"]) == 0.0)) {
        continue; // Skip electrons and non-charged ions
      }

      const Field3D Ni = getNoBoundary<Field3D>(species["density"]);
      const Field3D Ti = getNoBoundary<Field3D>(species["temperature"]);
      const BoutReal Mi = getNoBoundary<BoutReal>(species["AA"]);
      const BoutReal Zi = getNoBoundary<BoutReal>(species["charge"]);
      Field3D Vi = species.isSet("velocity")
        ? toFieldAligned(getNoBoundary<Field3D>(species["velocity"]))
        : zeroFrom(Ni);

      if (lower_y) {
        // Sum values, put result in mesh->ystart

        for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            auto i = indexAt(Ni, r.ind, mesh->ystart, jz);
            auto ip = i.yp();

            // Free gradient of log density and temperature
            // This ensures that the guard cell values remain positive
            // exp( 2*log(N[i]) - log(N[ip]) )

            const BoutReal Ni_im = limitFree(Ni[ip], Ni[i], density_boundary_mode);
            const BoutReal Ti_im = limitFree(Ti[ip], Ti[i], temperature_boundary_mode);

            // Calculate sheath values at half-way points (cell edge)
            const BoutReal nisheath = 0.5 * (Ni_im + Ni[i]);
            const BoutReal tisheath =
                floor(0.5 * (Ti_im + Ti[i]), 1e-5); // ion temperature

            // Sound speed squared
            BoutReal C_i_sq = (sheath_ion_polytropic * 2 * tisheath) / Mi;

            BoutReal visheath = -sqrt(C_i_sq);

            if (vD > 0){
              visheath *= (1+vD);
            }

            if (Vi[i] < visheath) {
              visheath = Vi[i];
            }

            ion_sum[i] -= Zi * nisheath * visheath;
          }
        }
      }

      if (upper_y) {
        // Sum values, put results in mesh->yend

        for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            auto i = indexAt(Ni, r.ind, mesh->yend, jz);
            auto im = i.ym();

            const BoutReal Ni_ip = limitFree(Ni[im], Ni[i], density_boundary_mode);
            const BoutReal Ti_ip = limitFree(Ti[im], Ti[i], temperature_boundary_mode);

            // Calculate sheath values at half-way points (cell edge)
            const BoutReal nisheath = 0.5 * (Ni_ip + Ni[i]);
            const BoutReal tisheath =
                floor(0.5 * (Ti_ip + Ti[i]), 1e-5); // ion temperature

            BoutReal C_i_sq = (sheath_ion_polytropic * 2 * tisheath) / Mi;

            BoutReal visheath = sqrt(C_i_sq);

            if (vD > 0){
              visheath *=(1-vD);
            }
            if (Vi[i] > visheath) {
              visheath = Vi[i];
            }

            ion_sum[i] += Zi * nisheath * visheath;
          }
        }
      }
    }

  //////////////////////////////////////////////////////////////////
  // Iterate through all ions
  for (auto& kv : allspecies.getChildren()) {
    if (kv.first == "e") {
      continue; // Skip electrons
    }

    Options& species = allspecies[kv.first]; // Note: Need non-const

    // Ion charge
    const BoutReal Zi = species.isSet("charge") ? get<BoutReal>(species["charge"]) : 0.0;

    if (Zi == 0.0) {
      continue; // Neutral -> skip
    }

    // Characteristics of this species
    const BoutReal Mi = get<BoutReal>(species["AA"]);

    // Density and temperature boundary conditions will be imposed (free)
    Field3D Ni = toFieldAligned(floor(getNoBoundary<Field3D>(species["density"]), 0.0));
    Field3D Ti = toFieldAligned(getNoBoundary<Field3D>(species["temperature"]));
    Field3D Pi = species.isSet("pressure")
      ? toFieldAligned(getNoBoundary<Field3D>(species["pressure"]))
      : Ni * Ti;

    // Get the velocity and momentum
    // These will be modified at the boundaries
    // and then put back into the state
    Field3D Vi = species.isSet("velocity")
      ? toFieldAligned(getNoBoundary<Field3D>(species["velocity"]))
      : zeroFrom(Ni);
    Field3D NVi = species.isSet("momentum")
      ? toFieldAligned(getNoBoundary<Field3D>(species["momentum"]))
      : Mi * Ni * Vi;

    if (lower_y) {
      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          auto i = indexAt(Ni, r.ind, mesh->ystart, jz);
          auto ip = i.yp();
          auto im = i.ym();

          // Free gradient of log electron density and temperature
          // This ensures that the guard cell values remain positive
          // exp( 2*log(N[i]) - log(N[ip]) )

          Ni[im] = limitFree(Ni[ip], Ni[i], density_boundary_mode);
          Ti[im] = limitFree(Ti[ip], Ti[i], temperature_boundary_mode);
          Pi[im] = limitFree(Pi[ip], Pi[i], pressure_boundary_mode);

          // Calculate sheath values at half-way points (cell edge)
          const BoutReal nisheath = 0.5 * (Ni[im] + Ni[i]);
          const BoutReal tisheath =
              floor(0.5 * (Ti[im] + Ti[i]), 1e-5); // ion temperature

          // Ion speed into sheath
          BoutReal C_i_sq = (sheath_ion_polytropic * 2 * tisheath) / Mi;

          BoutReal visheath = -sqrt(C_i_sq); // Negative -> into sheath

          if (vD > 0){
            visheath *= (1+vD);
          }

          if (Vi[i] < visheath) {
            visheath = Vi[i];
          }

          if (no_flow) {
            visheath = 0.0;
          }

          // Set boundary conditions on flows
          Vi[im] = 2. * visheath - Vi[i];
          NVi[im] = 2. * Mi * nisheath * visheath - NVi[i];
        }
      }
    }
    if (upper_y) {
      // Note: This is essentially the same as the lower boundary,
      // but with directions reversed e.g. ystart -> yend, ip <-> im
      //
      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          auto i = indexAt(Ni, r.ind, mesh->yend, jz);
          auto ip = i.yp();
          auto im = i.ym();

          // Free gradient of log electron density and temperature
          // This ensures that the guard cell values remain positive
          // exp( 2*log(N[i]) - log(N[ip]) )

          Ni[ip] = limitFree(Ni[im], Ni[i], density_boundary_mode);
          Ti[ip] = limitFree(Ti[im], Ti[i], temperature_boundary_mode);
          Pi[ip] = limitFree(Pi[im], Pi[i], pressure_boundary_mode);

          // Calculate sheath values at half-way points (cell edge)
          const BoutReal nisheath = 0.5 * (Ni[ip] + Ni[i]);
          const BoutReal tisheath =
              floor(0.5 * (Ti[ip] + Ti[i]), 1e-5); // ion temperature

          // Ion speed into sheath
          BoutReal C_i_sq = (sheath_ion_polytropic * 2 * tisheath) / Mi;

          BoutReal visheath = sqrt(C_i_sq); // Positive -> into sheath
          
          if (vD > 0){
              visheath *=(1-vD);
          }

          if (Vi[i] > visheath) {
            visheath = Vi[i];
          }

          if (no_flow) {
            visheath = 0.0;
          }

          // Set boundary conditions on flows
          Vi[ip] = 2. * visheath - Vi[i];
          NVi[ip] = 2. * Mi * nisheath * visheath - NVi[i];
        }
      }
    }
    // Finished boundary conditions for this species
    // Put the modified fields back into the state.

    Ni.clearParallelSlices();
    Ti.clearParallelSlices();
    Pi.clearParallelSlices();
    setBoundary(species["density"], fromFieldAligned(Ni));
    setBoundary(species["temperature"], fromFieldAligned(Ti));
    setBoundary(species["pressure"], fromFieldAligned(Pi));

    if (species.isSet("velocity")) {
      Vi.clearParallelSlices();
      setBoundary(species["velocity"], fromFieldAligned(Vi));
    }

    if (species.isSet("momentum")) {
      NVi.clearParallelSlices();
      setBoundary(species["momentum"], fromFieldAligned(NVi));
    }
  }
}
}

