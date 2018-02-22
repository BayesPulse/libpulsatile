#ifndef GUARD_ss_draw_sd_randomeffects_h
#define GUARD_ss_draw_sd_randomeffects_h

#include <RcppArmadillo.h>
#include <RInside.h>
#include <math.h>
#include "mh.h"
#include "patient.h"


//
// SS_DrawFixedEffects
//   Modified Metropolis Hastings sampler instantiating the mmh class for
//   sample the mean mass & width
//

class SS_DrawSDRandomEffects : public ModifiedMetropolisHastings<Patient, Patient, double, ProposalVariance>
{

  public:

    //
    // Constructor
    //   pass the proposal variance parameters to the constructor
    //
    SS_DrawSDRandomEffects(double in_pv, // double or arma::vec
                           int in_adjust_iter,
                           int in_max_iter,
                           double in_target_ratio) :
      ModifiedMetropolisHastings
      <Patient, Patient, double,
       ProposalVariance>::ModifiedMetropolisHastings(in_pv,
                                                     in_adjust_iter,
                                                     in_max_iter,
                                                     in_target_ratio) { };

  private:

    bool parameter_support(double val, Patient *patient) { 
        return (val > 0.0 && val < patient->priors->mass_sd_max);
    }

    double posterior_function(Patient *patient, double proposal, bool *notused) {


      return ;

    }

};

#endif

