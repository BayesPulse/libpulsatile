#ifndef GUARD_bpmod_singlesubject_draw_tvarscale_h
#define GUARD_bpmod_singlesubject_draw_tvarscale_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <math.h>
#include <bp_mcmc/mh.h>
#include <bp_datastructures/patient.h>


//
// SS_DrawTVarScale
//   Modified Metropolis Hastings sampler instantiating the mmh class for
//   sample the scale parameter for the t-distribution prior for drawing the
//   pulse mean/widths.
//

class SS_DrawTVarScale :
  public ModifiedMetropolisHastings<PulseEstimates, Patient, double, ProposalVariance>
{

  public:

    // Constructor
    SS_DrawTVarScale(double in_pv,
                     int in_adjust_iter,
                     int in_max_iter,
                     double in_target_ratio,
                     bool for_width,
                     bool verbose,
                     int verbose_iter) :
      ModifiedMetropolisHastings <PulseEstimates, Patient, double, ProposalVariance>::
      ModifiedMetropolisHastings(in_pv, in_adjust_iter, in_max_iter,
                                 in_target_ratio, verbose, verbose_iter) {

        // Choose which set of parameters to use: width or mass
        if (for_width) {
          est_mean_     = &PatientEstimates::width_mean;
          est_sd_       = &PatientEstimates::width_sd;
          randomeffect_ = &PulseEstimates::width;
          tvarscale_    = &PulseEstimates::tvarscale_width;
          parameter_name = "tvarscale width";
        } else {
          est_mean_     = &PatientEstimates::mass_mean;
          est_sd_       = &PatientEstimates::mass_sd;
          randomeffect_ = &PulseEstimates::mass;
          tvarscale_    = &PulseEstimates::tvarscale_mass;
          parameter_name = "tvarscale mass";
        }

      };

    // Pulse-specific estimate -- this function samples for each pulse
    void sample_pulses(Patient *patient, int iter) {

      for (auto &pulse : patient->pulses) {
        sample(&pulse, &(pulse.*tvarscale_), patient, iter);
      }

    }

  private:

    double PatientEstimates::*est_mean_;
    double PatientEstimates::*est_sd_;
    double PulseEstimates::*randomeffect_; //pulse specific mass or width
    double PulseEstimates::*tvarscale_;

    std::string parameter_name;
    std::string get_parameter_name() { return parameter_name; };

    bool parameter_support(double val, Patient *notused) {
      return (val > 0.0);
    };

    double posterior_function(PulseEstimates *pulse, double proposal, Patient *patient) {
      double old_gamma   = 0.0;
      double new_gamma   = 0.0;
      double prior_ratio = 0.0;
      double stdold      = 0.0;
      double stdnew      = 0.0;
      double re_old      = 0.0;
      double re_new      = 0.0;
      double re_ratio    = 0.0;
      PatientEstimates *est  = &patient->estimates;
      double patient_mean       = (*est).*est_mean_;
      double patient_sd         = (*est).*est_sd_;
      double pulse_randomeffect = (*pulse).*randomeffect_;
      double curr_scale         = (*pulse).*tvarscale_;

      // Shape, scale parameterized: 
      //    https://github.com/mmulvahill/r-source/blob/trunk/src/nmath/dgamma.c
      //    https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Distribution-functions
      old_gamma = Rf_dgamma(curr_scale, 2, 0.5, 0); 
      new_gamma = Rf_dgamma(proposal, 2, 0.5, 0);

      prior_ratio  = log(new_gamma) - log(old_gamma);

     // compute the normalization from the truncated t-distribution prior
      stdold       = (patient_mean / (patient_sd / sqrt(curr_scale));
      stdnew       = (patient_mean / (patient_sd / sqrt(proposal));
      re_old       = (pulse_randomeffect - patient_mean) * (pulse_randomeffect - patient_mean) * 0.5 * curr_scale;
      re_new       = (pulse_randomeffect - patient_mean) * (pulse_randomeffect - patient_mean)* 0.5 * proposal;
      re_ratio     = (re_old - re_new)/(patient_sd * patient_sd);
      re_ratio    += Rf_pnorm5(stdold, 0, 1, 1.0, 1.0) -  // second 1.0 does the log xform for us 
                     Rf_pnorm5(stdnew, 0, 1, 1.0, 1.0) -  // first 1.0 says to use lower tail      
                    log(curr_scale) + log(proposal); // the 1/2pi term in normal distirbution

      // Compute and acceptance ratio
      return prior_ratio + re_ratio;
    };

};




////
//// Defined functions for SD random effects MMH class
////
//
//// parameter_support()
////   Defines whether the proposal value is within the parameter support
//bool SS_DrawTVarScale::parameter_support(double val, Patient *notused) {
//  return (val > 0.0);
//}
//
//
//// posterior_function()
////   Calculates the acceptance ratio for use in modified metropolis hastings
////   sampler (inherited SS_DrawTVarScale::sample() function)
//double SS_DrawTVarScale::posterior_function(PulseEstimates *pulse, 
//                                            double proposal, 
//                                            Patient *patient) {
//
//  double old_gamma   = 0.0;
//  double new_gamma   = 0.0;
//  double prior_ratio = 0.0;
//  double stdold      = 0.0;
//  double stdnew      = 0.0;
//  double re_old      = 0.0;
//  double re_new      = 0.0;
//  double re_ratio    = 0.0;
//  PatientEstimates *est  = &patient->estimates;
//  double patient_mean       = (*est).*est_mean_;
//  double patient_sd         = (*est).*est_sd_;
//  double pulse_randomeffect = (*pulse).*randomeffect_;
//  double curr_scale         = (*pulse).*tvarscale_;
//
//  // Shape, scale parameterized: 
//  //    https://github.com/mmulvahill/r-source/blob/trunk/src/nmath/dgamma.c
//  //    https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Distribution-functions
//  old_gamma = Rf_dgamma(curr_scale, 2, 0.5, 0); 
//  new_gamma = Rf_dgamma(proposal, 2, 0.5, 0);
//
//  prior_ratio  = log(new_gamma) - log(old_gamma);
//
//  stdold       = (pulse_randomeffect) / (patient_sd / sqrt(curr_scale));
//  stdnew       = (pulse_randomeffect) / (patient_sd / sqrt(proposal));
//  re_old       = (pulse_randomeffect - patient_mean) * 0.5 * re_old * curr_scale;
//  re_new       = (pulse_randomeffect - patient_mean) * 0.5 * re_new * proposal;
//  re_ratio     = (re_old - re_new) / (patient_sd * patient_sd);
//  re_ratio    += Rf_pnorm5(stdold, 0, 1, 1.0, 1.0) -  // second 1.0 does the log xform for us 
//                 Rf_pnorm5(stdnew, 0, 1, 1.0, 1.0) -  // first 1.0 says to use lower tail      
//                 0.5 * log(curr_scale) + 0.5 * log(proposal); // the 1/2pi term in normal distirbution
//
//  // Compute and acceptance ratio
//  return prior_ratio + re_ratio;
//
//};

#endif


