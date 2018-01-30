#ifndef GUARD_ss_drawlocations_h
#define GUARD_ss_drawlocations_h

#include <RcppArmadillo.h>
#include <RInside.h>
#include <math.h>
#include "mh.h"
#include "patient.h"


// 
// SS_DrawLocations
//   Modified Metropolis Hastings sampler instantiating the mmh class for
//   sampling the pulse locations
//

class SS_DrawLocations : public ModifiedMetropolisHastings<Patient, double, ProposalVariance>
{

  public:
    // Constructors

  private:

    bool parameter_support(double val) { return true; }

    double posterior_function(Patient *patient, double proposal) {



      return 0.0; 

    }

}


#endif

