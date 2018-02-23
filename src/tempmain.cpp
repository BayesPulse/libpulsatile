
#include <RcppArmadillo.h>
#include <RInside.h>
//#include "catch.h"
//#include "counter.h"
//#include "proposalvariance.h"
#include "mh.h"
#include "patient.h"
#include "population.h"
#include "utils.h"
#include "ss_draw_fixedeffects.h"
#include "ss_draw_sdrandomeffects.h"
#include "ss_draw_baselinehalflife.h"
#include "ss_draw_locations.h"
#include "ss_draw_randomeffects.h"
#include "ss_draw_tvarscale.h"

int main(int argc, char **argv) {

  std::cout << "Hello world\n";

  // Create R instance (for RNGs) and set seed for reproducing test results
  RInside R;
  PulseUtils pu;
  pu.set_seed(171227);

  //// Create objects for checking/testing RNG implementation
  //arma::vec initial_means = { 2, 3 };
  //arma::vec initial_pvs   = { 0.7, 0.1 };

  //ProposalVariance2p pv(initial_pvs, 500, 25000, 0.25);
  //arma::mat cholmat = pv.getpsd();
  //double cholcov    = cholmat(0, 1);

  //// Test/Print RNG results
  //std::cout << "PV = " << pv.getpv() << "\n";
  //std::cout << "Choldecomp of PV (i.e. PSD, upper triangle form)  = " << pv.getpsd() << "\n";
  //std::cout << "PSD test that actually upper triangular form = " << cholcov << "\n";
  //std::cout << "rmvnorm = " << pu.rmvnorm(initial_means, pv.getpsd()) << "\n";

  ////double a = 11;
  ////arma::vec testvec = initial_means;
  ////testvec += a;
  ////std::cout << "vector += scalar: " << initial_means << " += " << a <<  " = " << testvec << "\n";
  ////// TRUE

  //// Testing using R objects in c++
  //NumericVector time(100); //, conc(100);
  //NumericVector conc     = rnorm(100, 3, 0.1);
  //NumericVector response = conc + rnorm(100, 1, 0.001);
  //for (int i = 0; i < time.size(); i++) time(i) = i + 1;


  //std::cout << "Time meta: size=" << time.size() << ", begin=" << time(0) << ", end=" << time(time.size() - 1) << "\n";

  ////std::cout << "NumericVector w/ values from 1 to 100 and 100 random normals and 100 more\n";
  ////for (int i = 0; i < time.size(); i++) {
  ////  std::cout << time(i) << ", " << conc(i) << ", " << response(i) << "\n";
  ////}

  //////////////// BEGIN LOADING DATA STRUCTURES ////////////////
  //               Data and Patient object setup               //
  ///////////////////////////////////////////////////////////////

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

  // Create patient data object
  PatientData pdone(thistime, conc);
  //Create priors object
  PatientPriors ppsingle(1.5, 100, 45, 100, 3.5, 100, 30, 100,
                         10, 100, 1000, 1000, 12, 0, 40);
  // Create estimates object (w/ starting vals)
  PatientEstimates pesingle(2.6, 45, 0.05, 3.5, 30, 12, 10, 10);

  // Create pointers
  PatientData * data = &pdone;
  PatientPriors * priors = &ppsingle;
  PatientEstimates * estimates = &pesingle;

  // Now take all of this and create a Patient object
  Patient pat(data, priors, estimates);
  Patient * patient = &pat;

  // Now add pulse estimates
  // Delete initial pulse for these tests
  ++pat.piter;
  pat.piter = pat.pulses.erase(pat.piter);

  // Now add a bunch of pulses (from matt's sim data)
  arma::vec location        { -37.65204, 41.04917, 125.58297, 306.32966,
                              461.65469, 616.95867, 835.64466, 1000.29703,
                              1149.08360, 1319.19888, 1414.31830 };
  arma::vec mass            { 6.635627,  5.940251, 13.463474,  7.649251,
                              4.309497, 12.285200,  4.241401,  3.734942,
                              3.553504,  4.589623, 7.364310 };
  arma::vec width           { 4.596894, 6.740636, 5.986173, 4.476816, 4.701983,
                              3.829487, 5.917652, 4.932739, 6.367862, 4.915543,
                              5.692852 };
  arma::vec tvarscale_mass  { 1.3362258, 0.9683636, 0.2660084, 0.1090139,
                              0.5613391, 0.4339317, 0.3184480, 0.2060616,
                              1.4164371, 2.1759781, 1.0744991 };
  arma::vec tvarscale_width { 1.0301410, 0.7512098, 1.9517880, 1.9834742,
                              1.2743407, 0.4760933, 1.7499186, 0.6057495,
                              1.5219036, 0.8946568, 0.8632652 };

  for (int i = 0; i < location.n_elem; i++) {
    PulseEstimate pulse(location(i), mass(i), width(i), tvarscale_mass(i), tvarscale_width(i),
                        pat.estimates->get_decay(), pat.data->time);
    pat.pulses.push_back(pulse);
  }
  ////////////////////////////////////////////////////////////
  //////////////// END LOADING DATA STRUCTURES ///////////////
  ////////////////////////////////////////////////////////////


  // Create sampler object 
  arma::vec bhl_pv = { 0.5, 45 };
  SS_DrawFixedEffects draw_fixed_effects(1.1, 500, 25000, 0.35);
  SS_DrawSDRandomEffects draw_sd_pulse_masses(2, 500*11, 25000*11, 0.35);
  SS_DrawBaselineHalflife draw_baselinehalflife(bhl_pv, 500, 25000, 0.25);
  SS_DrawLocationsStrauss draw_pulse_locations_strauss(10, 500*11, 25000*11, 0.35);
  SS_DrawRandomEffects draw_pulse_masses(1.1, 500*11, 25000*11, 0.35);
  SS_DrawTVarScale draw_pulse_tvarscale(1.01, 500*11, 25000*11, 0.35);


  // Sample MMH objects
  for (int i = 0; i < 250000; i++) {
    draw_fixed_effects.sample(patient, &patient->estimates->mass_mean);
    draw_sd_pulse_masses.sample(patient, &patient->estimates->mass_sd, patient);
    draw_baselinehalflife.sample(patient, &patient->estimates->baseline_halflife);
    draw_pulse_locations_strauss.sample_pulses(patient);
    draw_pulse_masses.sample_pulses(patient);
    draw_pulse_tvarscale.sample_pulses(patient);
  }

  return 0;

}


