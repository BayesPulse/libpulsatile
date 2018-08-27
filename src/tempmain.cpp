#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RInside.h>
#include <singlesubject.h>
#include <bp_mcmc/utils.h>

using namespace Rcpp;

int main(int argc, char **argv) {

  std::cout << "Hello world\n";

  // Create R instance (for RNGs) and set seed for reproducing test results
  RInside R;
  PulseUtils pu;
  pu.set_seed(999999);


  //////////////// BEGIN LOADING DATA STRUCTURES ////////////////

  // Create patient data
  Rcpp::NumericVector thistime(144);
  for (int i = 0; i < thistime.size(); i++)  thistime(i) = (i + 1) * 10;
  Rcpp::NumericVector conc =
  { 5.237937,  5.156800,  5.773434,  6.364134,  9.872640,  8.782145, 8.633523,
    7.166886,  6.063834,  6.393283,  5.547755,  5.609554, 17.809125,
    14.864233, 14.488007, 11.012978, 10.074941,  8.579306,  7.925953,
    7.505312,  7.408293, 6.296038,  6.138372,  6.114483,  4.589039,  4.787283,
    4.125096,  4.560712, 3.838159,  3.611227, 10.149188,  9.114352,  8.475730,
    8.166254,  7.115865, 6.364551,  6.076659,  4.936654,  4.815760,  5.000790,
    4.405922,  4.621510, 3.612532,  3.861247,  3.480776,  4.382934,  6.908732,
    5.677212,  5.697563, 5.049410,  5.033481,  5.036818,  4.766425,  4.108461,
    4.493277,  3.970064, 3.479031,  3.105488,  3.182802,  2.832481,  2.904569,
    13.393956, 12.160363, 12.604881,  9.282701,  8.971699,  7.872894,
    7.341750,  6.129950,  6.754877, 6.236324,  5.221481,  4.712401,  4.251457,
    4.107650,  3.642404,  3.409108, 3.673586,  3.573048,  3.006202,  3.574850,
    2.984541,  2.942794,  6.142817, 6.239127,  5.607103,  5.781841,  5.365100,
    4.483500,  4.318313,  3.729648, 4.186974,  3.661007,  3.587046,  3.420915,
    3.554909,  2.951518,  3.197826, 2.768718,  4.050486,  6.375893,  5.984611,
    4.975292,  4.758750,  5.006260, 4.137620,  3.676703,  4.133922,  3.445431,
    3.813060,  3.304948,  3.308399, 3.362024,  3.239059,  5.554801,  5.157033,
    5.049831,  5.366527,  4.410822, 4.798887,  4.386467,  4.013058,  3.749040,
    3.580090,  3.434933,  3.086843, 3.471470,  3.334008,  2.913723,  2.980340,
    2.862694,  5.891940,  6.136203, 6.251388,  6.170612,  5.191264,  4.810252,
    4.989489,  3.951583,  3.793401, 4.320815, 10.128100,  9.961510,  8.155648 };

  //
  // Create priors, starting values, and proposal variances objects
  //
  Rcpp::List priors = List::create(Named("baseline_mean")           = 2.6,
                                   Named("baseline_variance")       = 100,
                                   Named("halflife_mean")           = 45,
                                   Named("halflife_variance")       = 100,
                                   Named("mass_mean")               = 3.5,
                                   Named("mass_variance")           = 100,
                                   Named("width_mean")              = 30,
                                   Named("width_variance")          = 100,
                                   Named("mass_sdmax")            = 1000,
                                   Named("width_sdmax")           = 1000,
                                   Named("error_alpha")             = 0.0001,
                                   Named("error_beta")              = 0.0001,
                                   Named("pulse_count")             = 12,
                                   Named("location_prior_type")     = "strauss",
                                   Named("strauss_repulsion")       = 0,
                                   Named("strauss_repulsion_range") = 40);
  priors.attr("class") = "bp_priors";

  // Create estimates object (w/ starting vals)
  Rcpp::List startingvals = List::create(Named("baseline")    = 2.6,
                                         Named("halflife")    = 45,
                                         Named("errorsq")     = 0.5,
                                         Named("mass_mean")   = 3.5,
                                         Named("width_mean")  = 30,
                                         Named("pulse_count") = 12,
                                         Named("mass_sd")     = 10,
                                         Named("width_sd")    = 10);
  startingvals.attr("class") = "bp_startingvals";

  Rcpp::List proposalvars = List::create(Named("mass_mean") = 1.1,
                                         Named("width_mean") = 1.1,
                                         Named("mass_sd") = 2,
                                         Named("width_sd") = 2,
                                         Named("baseline") = 0.5,
                                         Named("halflife") = 45,
                                         Named("location") = 1000,
                                         Named("pulse_mass") = 1,
                                         Named("pulse_width") = 10 ,
                                         Named("sdscale_pulse_mass") = 1,
                                         Named("sdscale_pulse_width") = 1);
  proposalvars.attr("class") = "bp_proposalvariance";

  //////////////// END LOADING DATA STRUCTURES ///////////////


  // Create sampler object 
  Rcpp::List rtn_list;
  rtn_list = singlesubject_(conc, thistime, priors, proposalvars, startingvals,
                            //10000, 50, 10000, true, 500, 25000, 0.25, 0.35);
                            100, 1, 1, true, 500, 25000, 0.25, 0.35);

  return 0;

}


