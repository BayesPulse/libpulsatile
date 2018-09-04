#ifndef GUARD_bpmod_singlesubject_draw_locations_h
#define GUARD_bpmod_singlesubject_draw_locations_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <math.h>
#include <bp_datastructures/patient.h>
#include <bp_mcmc/mh.h>
#include <bp_mcmc/utils.h>


// 
// SS_DrawLocation
//   Base class for pulse locations -- 
//   Modified Metropolis Hastings sampler instantiating the mmh class for
//   sampling the pulse locations
//

class SS_DrawLocations : public ModifiedMetropolisHastings<PulseIter, Patient, double, ProposalVariance>
{

  public:
    // Constructors
    SS_DrawLocations(double in_pv, int in_adjust_iter, int in_max_iter,
                     double in_target_ratio, bool in_verbose, int in_verbose_iter) :
      ModifiedMetropolisHastings <PulseIter, Patient, double, ProposalVariance> :: 
      ModifiedMetropolisHastings(in_pv, in_adjust_iter, in_max_iter,
                                 in_target_ratio, in_verbose, in_verbose_iter) { };

    virtual ~SS_DrawLocations() { }

    std::string parameter_name = "pulse locations";
    std::string get_parameter_name() { return parameter_name; };

    // Pulse level estimates need to be done at the pulse level
    void sample_pulses(Patient *patient, int iter) {

      for (auto pulse = patient->pulses.begin(); pulse != patient->pulses.end(); ++pulse) {
        sample(&pulse, &(pulse->time), patient, iter);
      }

    }

};


#endif

