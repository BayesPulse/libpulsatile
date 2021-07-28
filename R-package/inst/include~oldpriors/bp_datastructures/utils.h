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
        { 2.6814358424853, 3.61930378208367, 4.5835959631532, 5.39164025413596,
          4.60258017744968, 5.83757715322813, 5.09375858456153, 4.46036157711056,
          4.22856819090718, 3.74957838153654, 3.75609077429634, 3.79411659707651,
          3.45805260092717, 3.61457289001595, 2.97937379998757, 3.23557654740535,
          3.37949836675526, 3.08395676488862, 3.29942661590539, 3.467284650064,
          2.96904974806636, 3.4012276452592, 3.14614619252082, 2.93431531436762,
          3.04134596753656, 3.16108792824359, 2.74099223041317, 2.8699720051147,
          3.23658024143572, 5.19212131762849, 5.75254811242479, 5.49677909914257,
          5.52635482285921, 4.94530043569383, 4.93650498103998, 7.31260080328799,
          7.12839518430793, 6.80174472098413, 7.63568878301611, 6.24612828790289,
          6.03269015404348, 4.70255921360491, 5.03914919538947, 4.21976716592007,
          4.5274923600401, 4.0901150408045, 3.45540128735323, 3.18741144528796,
          4.27177037249524, 6.863361313599, 5.95369395309416, 6.45675517565742,
          6.77894223698634, 4.98197746165615, 4.83024785565837, 4.33655737698785,
          4.2962112048339, 4.23747662113222, 3.42746534312757, 3.66541920773442,
          2.97883689544759, 3.40956869822447, 2.64176227029741, 3.6897393626309,
          3.18030179778705, 3.49793997046746, 3.0979109277295, 3.20592508831985,
          4.96256037135427, 5.55992818548297, 6.37703530176968, 4.97415141089928,
          5.09490646872372, 4.86449992541152, 4.85230880348274, 5.3152210319564,
          5.49178722081631, 5.37978623379431, 4.84506978201334, 4.83524164889408,
          4.5230894407225, 4.21198537611242, 4.74196595628771, 6.41071424973019,
          7.96640733345621, 7.23367560511177, 6.29354138548944, 5.80755330229538,
          5.62640845515504, 4.68537869434213, 4.97610406855091, 4.92376102540039,
          5.61631365195449, 4.95454646496695, 4.31629638822296, 4.4493813988051,
          4.03561211218328, 5.93303717625482, 5.46421418220275, 5.14575060734008,
          5.20019144472861, 4.55307603773356, 4.42996713679056, 3.91582994527552,
          3.96257479535214, 3.41896505108181, 3.3348628406967, 3.17450031947394,
          3.4093277142779, 2.82261506453114, 3.29827697890822, 2.42123259143846,
          3.41368318065368, 2.85054712083439, 3.11556223300741, 2.71361603688389,
          2.94197971715263, 2.88786583119099, 2.98076624150669, 3.62782435408637,
          4.62560519056033, 4.46845075720742, 4.81575824835943, 3.98543574808699,
          3.46347068453034, 3.68228608222232, 3.53695782519207, 3.56394158512654,
          3.55280962733125, 3.75170861906579, 2.93317036381674, 4.23415770057993,
          4.71644511719147, 4.04372677190788, 4.32006371035954, 3.97229929496866,
          3.85122465873892, 4.05022084815804, 3.19514276239586, 3.16839853737732,
          3.0116536403404, 2.72138430008304, 3.27921071119721, 3.079000425982 };

      // Create patient data object
      PatientData data(thistime, conc);
      //Create priors object
      PatientPriors priors(2.6, 100, 45, 100, 3.5, 100, 42, 1000,
                           100, 1000, 0.0001, 0.0001, 12, 0, 40);
      // Create estimates object (w/ starting vals)
      PatientEstimates estimates(2.6, 45, 0.005, 3.5, 42, 1.6, 35);

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
      arma::vec location        { 26.5415235605463, 174.639932001226,
        298.621168723401, 360.553287402341, 494.611547829321, 689.092415108905,
        763.890165489543, 839.800270423608, 925.802510831759, 975.473197559335,
        1199.00866308622, 1322.82471326204 };
      arma::vec mass            { 4.24073843802133, 0.765557146882449,
        4.08258116960848, 4.76556854871956, 4.79030238300445, 4.08353146779168,
        2.35870955639405, 5.37561799380853, 1.26164550948196, 2.73323761200587,
        2.43613946981746, 2.19844374502296 };
      arma::vec width           { 88.1111213719339, 110.366726051791,
        36.036933446449, 49.107878599033, 20.2837571325304, 36.4618443022452,
        57.3848450667396, 57.4019670930382, 1.92996218540478, 15.4293652120993,
        139.865750584882, 43.9611972363173 };
      arma::vec tvarscale_mass  { 0.848675552152578, 0.710826927142778,
        0.904553732882957, 0.420310565275885, 0.841377131089954,
        0.456265805712836, 1.2856220853906, 1.23367393752701, 0.672893954384988,
        0.940914640905448, 1.60380806975763, 0.631232048387259 };
      arma::vec tvarscale_width { 0.680352135467811, 0.616970113498595,
        0.748687310834558, 0.465599040614681, 1.33113209570755,
        0.548808523542267, 0.62902939160845, 1.42678455374602,
        0.669347263951419, 0.576907900295665, 0.587171571105131,
        0.85079857539479 };

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