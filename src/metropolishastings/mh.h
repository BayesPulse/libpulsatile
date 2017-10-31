// Abstract class for Metropolis Hastings samplers
// Demonstrates declaration of a constructors and destructor for the Cat class

//#include <iostream.h> // for cout
//#include <Rcpp.h>
//using namespace Rcpp;


template <class T_myMMH>
class ModifiedMetropolisHastings
{

  public:
    // sample from posterior
    //   - runs 1 iteration
    //   - inputs/outputs or changes internally?
    double sample() {
      return static_cast<T_myMMH &>(*this).sample();
      // or
      //return static_cast<T_myMMH *>(this)->sample();
    };

  private:
    double draw_proposal(double current, double pv.getpv());
    double parameter_support();  // i.e. truncation logic
    double posterior_function(); // logrho_calculation
    Proposalvariance pv;          // this should be an object of this class
};


double ModifiedMetropolisHastings::draw_proposal(double current, double proposal_sd)
{
  return Rf_rnorm(current, proposal_sd);
}


















class sd_widths : ModifiedMetropolisHastings<sd_widths> {

};


class sd_mass : ModifiedMetropolisHastings<sd_widths> {

};

class baseline : ModifiedMetropolisHastings<sd_widths> {

};

class halflife : ModifiedMetropolisHastings<sd_widths> {

};

template <class T> 
double sample(base<T> & mmh) {

  double proposal = mmh.draw_proposal;
  if (mmh.parameter_support()) {


  } else {
    return
  }


}





template <class T_myMMH>
class ModifiedMetropolisHastings 
{

  public: 
    double
  // template method
  //  how to choose return type in derived class?
  myMMH draw_proposal() {
    static_cast<myMMH *>(this)->draw_proposal();
  };
  myMMH within_parameter_support() {
    static_cast<myMMH *>(this)->within_parameter_support();
  };
  myMMH within_parameter_support() {
    static_cast<myMMH *>(this)->within_parameter_support();
  };

};




