#ifndef GUARD_bp_datastructures_utils_h
#define GUARD_bp_datastructures_utils_h

#include <RcppArmadillo.h>
#include <bp_datastructures/patient.h>
#ifndef NORINSIDE
#include <RInside.h>
#endif


//
// utils.h
//    Defines the DataStructuresUtils class, containing miscellaneous functions
//    and calculations needed for testing and using the bp datastructures.
//

class DataStructuresUtils {

  public:

    // create_new_test_patient_obj()
    //   Function for initializing a patient object with a default dataset and
    //   specs
    //
    Patient create_new_test_patient_obj()
    {

      // Create patient data
      Rcpp::NumericVector thistime(144);
      for (int i = 0; i < thistime.size(); i++)  thistime(i) = (i + 1) * 10;
      Rcpp::NumericVector conc =
      { 2.681436, 3.619304, 4.583596, 5.391640, 4.602580, 5.837577, 5.093759,
        4.460362, 4.228568, 3.749578, 3.756091, 3.794117, 3.458053, 3.614573,
        2.979374, 3.235577, 3.379498, 3.083957, 3.299427, 3.467285, 2.969050,
        3.401228, 3.146146, 2.934315, 3.041346, 3.161088, 2.740992, 2.869972,
        3.236580, 5.192121, 5.752548, 5.496779, 5.526355, 4.945300, 4.936505,
        7.312601, 7.128395, 6.801745, 7.635689, 6.246128, 6.032690, 4.702559,
        5.039149, 4.219767, 4.527492, 4.090115, 3.455401, 3.187411, 4.271770,
        6.863361, 5.953694, 6.456755, 6.778942, 4.981977, 4.830248, 4.336557,
        4.296211, 4.237477, 3.427465, 3.665419, 2.978837, 3.409569, 2.641762,
        3.689739, 3.180302, 3.497940, 3.097911, 3.205925, 4.962560, 5.559928,
        6.377035, 4.974151, 5.094906, 4.864500, 4.852309, 5.315221, 5.491787,
        5.379786, 4.845070, 4.835242, 4.523089, 4.211985, 4.741966, 6.410714,
        7.966407, 7.233676, 6.293541, 5.807553, 5.626408, 4.685379, 4.976104,
        4.923761, 5.616314, 4.954546, 4.316296, 4.449381, 4.035612, 5.933037,
        5.464214, 5.145751, 5.200191, 4.553076, 4.429967, 3.915830, 3.962575,
        3.418965, 3.334863, 3.174500, 3.409328, 2.822615, 3.298277, 2.421233,
        3.413683, 2.850547, 3.115562, 2.713616, 2.941980, 2.887866, 2.980766,
        3.627824, 4.625605, 4.468451, 4.815758, 3.985436, 3.463471, 3.682286,
        3.536958, 3.563942, 3.552810, 3.751709, 2.933170, 4.234158, 4.716445,
        4.043727, 4.320064, 3.972299, 3.851225, 4.050221, 3.195143, 3.168399,
        3.011654, 2.721384, 3.279211, 3.079000 };




      // Create patient data object
      PatientData data(thistime, conc);
      //Create priors object
      PatientPriors priors(2.6, 100, 45, 100, 3.5, 100, 42, 100,
                           10, 150, 1000, 1000, 12, 0, 40);
      // Create estimates object (w/ starting vals)
      PatientEstimates estimates(2.6, 45, 0.05, 3.5, 30, 10, 10);

      // Now take all of this and create a Patient object
      Patient pat(data, priors, estimates);

      return pat;

    }


    Patient* add_default_pulses(Patient * patient) {

      // Now add pulse estimates
      // Delete initial pulse for these tests
      ++patient->piter;
      patient->piter = patient->pulses.erase(patient->piter);

      // Now add a bunch of pulses (from matt's sim data)
      arma::vec location        { 26.54152, 174.63993, 298.62117, 360.55329,
        494.61155, 689.09242, 763.89017, 839.80027, 925.80251, 975.47320,
        1199.00866, 1322.82471 };
      arma::vec mass            { 4.2407384, 0.7655571, 4.0825812, 4.7655685,
        4.7903024, 4.0835315, 2.3587096, 5.3756180, 1.2616455, 2.7332376,
        2.4361395, 2.1984437 };
      arma::vec width           { 88.111121, 110.366726, 36.036933, 49.107879,
        20.283757, 36.461844, 57.384845, 57.401967, 1.929962, 15.429365,
        139.865751, 43.961197 };
      arma::vec tvarscale_mass  { 0.8486756, 0.7108269, 0.9045537, 0.4203106,
        0.8413771, 0.4562658, 1.2856221, 1.2336739, 0.6728940, 0.9409146,
        1.6038081, 0.6312320 };
      arma::vec tvarscale_width { 0.6803521, 0.6169701, 0.7486873, 0.4655990,
        1.3311321, 0.5488085, 0.6290294, 1.4267846, 0.6693473, 0.5769079,
        0.5871716, 0.8507986 };

      for (arma::uword i = 0; i < location.n_elem; i++) {
        PulseEstimates pulse(location(i), mass(i), width(i), tvarscale_mass(i), tvarscale_width(i),
                             patient->estimates.get_decay(), patient->data.time);
        patient->pulses.push_back(pulse);
      }

      return patient;

    }

};

//----------------------------------------------------------------------------//
// END OF FILE
//----------------------------------------------------------------------------//

#endif
