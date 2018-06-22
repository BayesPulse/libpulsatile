#include <RcppArmadillo.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif
#include <bp_datastructures/chains.h>
#include <testing/catch.h>


//
// chains_tests.cpp
//     Test chain class & functions
//

//

TEST_CASE( "first chains test -- single subject, single hormone", "[chains]" ) {

  ////////////////// BEGIN  LOADING DATA STRUCTURES ///////////////
  ////               Data and Patient object setup               // 
  /////////////////////////////////////////////////////////////////

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
  PatientEstimates pesingle(2.6, 45, 0.05, 3.5, 30, 10, 10);

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
    PulseEstimates pulse(location(i), mass(i), width(i), tvarscale_mass(i), tvarscale_width(i),
                         pat.estimates->get_decay(), pat.data->time);
    pat.pulses.push_back(pulse);
  }
  ////////////////// END LOADING DATA STRUCTURES ///////////////


  // Single subject constructor
  Chains chains(100, 2, 10, false);

  SECTION( "Check constructor creates empty chain objects" ) {
    REQUIRE( chains.patient_chain.n_rows == 45 );
    REQUIRE( chains.patient_chain.n_cols == 9 );
    REQUIRE( std::abs(chains.patient_chain.max() - chains.patient_chain.min()) < 0.000001 );
    REQUIRE( chains.patient_chain(11, 1) == 0.0 );
    REQUIRE( chains.patient_chain(12, 1) == 0.0 );
  }

  chains.save_sample(&pat, 30);
  chains.save_sample(&pat, 32);
  chains.save_sample(&pat, 52);

  SECTION( "Check save_sample() function populates common chain" ) {
    REQUIRE( std::abs(chains.patient_chain.max() - chains.patient_chain.min()) > 0.000001 );
    REQUIRE( chains.patient_chain(10, 0) == 30 );
    REQUIRE( chains.patient_chain(10, 1) == 11 );
    REQUIRE( chains.patient_chain(10, 2) == 2.6 );
    REQUIRE( chains.patient_chain(10, 3) == 3.5 );
    REQUIRE( chains.patient_chain(10, 4) == 30 );
    REQUIRE( chains.patient_chain(10, 5) == 45 );
    REQUIRE( chains.patient_chain(10, 6) == 0.05 );
    REQUIRE( chains.patient_chain(10, 7) == 10 );
    REQUIRE( chains.patient_chain(10, 8) == 10 );
    REQUIRE( chains.patient_chain(11, 0) == 32 );
    REQUIRE( chains.patient_chain(11, 1) == 11 );
    REQUIRE( chains.patient_chain(21, 0) == 52 );
    REQUIRE( chains.patient_chain(21, 1) == 11 );
  }

  //std::cout << "chain 2 is not empty " << chains.pulse_chains[2] << std::endl;
  //std::cout << "chain 3 is empty " << (chains.pulse_chains[3] == chains.pulse_chains.end()) << std::endl;
  SECTION( "Check save_sample() function populates pulse chain vector of matrices" ) {
    REQUIRE( chains.pulse_chains[0].max() - chains.pulse_chains[0].min() > 1000 );
    REQUIRE( chains.pulse_chains[1].max() - chains.pulse_chains[1].min() > 1000);
    REQUIRE( chains.pulse_chains[2].max() - chains.pulse_chains[2].min() > 1000 );
    REQUIRE( &chains.pulse_chains[2] == &chains.pulse_chains.back() );
  }



}
