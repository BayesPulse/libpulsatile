#ifndef GUARD_patientestimates_h
#define GUARD_patientestimates_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <bp_mcmc/utils.h>

//
// populationestimates.h
//   defining the estimates object used at the population level
//
// Author: Max McGrath
// Notes:
//

using namespace Rcpp;

struct PopulationEstimates {
  double mass_mean;
  double mass_s2s_sd;
  double mass_p2p_sd;

  double width_mean;
  double width_s2s_sd;
  double width_p2p_sd;

  double baseline_mean;
  double baseline_sd;

  double halflife_mean;
  double halflife_sd;

  PopulationEstimates(double sv_mass_mean,
                      double sv_mass_s2s_sd,
                      double sv_mass_p2p_sd,
                      double sv_width_mean,
                      double sv_width_s2s_sd,
                      double sv_width_p2p_sd,
                      double sv_baseline_mean,
                      double sv_baseline_sd,
                      double sv_halflife_mean,
                      double sv_halflife_sd) {

  mass_mean = sv_mass_mean;
  mass_s2s_sd = sv_mass_s2s_sd;
  mass_p2p_sd = sv_mass_p2p_sd;
  width_mean = sv_width_mean;
  width_s2s_sd = sv_width_s2s_sd;
  width_p2p_sd = sv_width_p2p_sd;
  baseline_mean = sv_baseline_mean;
  baseline_sd = sv_baseline_sd;
  halflife_mean = sv_halflife_mean;
  halflife_sd = sv_halflife_sd;

  }


};


#endif
