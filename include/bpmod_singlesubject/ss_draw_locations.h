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
    SS_DrawLocations(double in_pv, // double or arma::vec
                     int in_adjust_iter,
                     int in_max_iter,
                     double in_target_ratio) :
      ModifiedMetropolisHastings
      <PulseIter, Patient, double,
       ProposalVariance>::ModifiedMetropolisHastings(in_pv,
                                                     in_adjust_iter,
                                                     in_max_iter,
                                                     in_target_ratio) { };
    virtual ~SS_DrawLocations() { }
    // Pulse level estimates need to be done at the pulse level
    void sample_pulses(Patient *patient, int iter) {

      //PulseIter pulse = patient->pulses.begin();
      //PulseConstIter pulse_end = patient->pulses.end();

      //while (pulse != pulse_end) {
      //  // Sample pulse,
      //  //   note: &(*pulse) derefs iter, then gets address of underlying obj
      //  sample(&(*pulse), &pulse->time, patient, iter);
      //  pulse++;
      //}
      //for (auto &pulse : patient->pulses) {
      for (auto pulse = patient->pulses.begin(); pulse != patient->pulses.end(); ++pulse) {
        //sample(&(*pulse), &pulse->time, patient, iter);
        sample(&pulse, &(pulse->time), patient, iter);
      }

    }

};


#endif

