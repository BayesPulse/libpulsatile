#ifndef GUARD_bpmod_singlesubject_draw_error_h
#define GUARD_bpmod_singlesubject_draw_error_h

#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <math.h>
#include <bp_datastructures/patient.h>

//
// SS_DrawError
//   Draw the model error variance from the inverse Gamma distribution with Gibbs
//   Modeling variance not precision.
//

class SS_DrawError 
{

  public:
    // Constructors
    SS_DrawError() { };

    void sample(Patient *patient) {

      int N        = patient->data.number_of_obs;
      double alpha = patient->priors.error_alpha;
      double beta  = patient->priors.error_beta;
      double ssq   = patient->get_sumerrorsquared(false);
      //Rcpp::Rcout << "SSQ:" << ssq << std::endl;

      //patient->estimates.errorsq  = 1 / Rf_rgamma(alpha + N / 2, beta + (0.5 * ssq));
      patient->estimates.errorsq  = 1 / Rf_rgamma(alpha + N / 2, 1 / (1 / beta + 0.5 * ssq));
      //parms->sigma  = 1 / Rf_rgamma(priors->err_alpha + N / 2, 1 / (1 / priors->err_beta + 0.5 * ssq));
    }

};

#endif

