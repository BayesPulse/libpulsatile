// Abstract class for Metropolis Hastings samplers
// Demonstrates declaration of a constructors and destructor for the Cat class

//#include <iostream.h> // for cout
//#include <Rcpp.h>
//using namespace Rcpp;


template <class T> 
class ModifiedMetropolisHastings
{

  public: 
    // constructor and destructor
    ModifiedMetropolisHastings(double starting_value,
                               double proposal_stddev);
    ~ModifiedMetropolisHastings();
    // public methods
    double sample(); // sample from posterior
                     //   - runs 1 iteration
                     //   - inputs/outputs or changes internally?

  private:
    double posterior_function(); // logrho_calculation
    double parameter_support();  // i.e. truncation logic
    Proposalvariance pv;          // this should be an object of this class

    Posterior posterior_function;
    Support parameter_support;
    MCMC_state current_state;

};


class ModifiedMetropolisHastings::










// constructor of Cat,
Cat::Cat(int initialAge)
{
  itsAge = initialAge;
}

Cat::~Cat()                 // destructor, takes no action
{
}

// GetAge, Public accessor function
// returns value of itsAge member
int Cat::GetAge()
{
  return itsAge;
}

// Definition of SetAge, public
// accessor function

void Cat::SetAge(int age)
{
  // set member variable its age to
  // value passed in by parameter age
  itsAge = age;
}

// definition of Meow method
// returns: void
// parameters: None
// action: Prints "meow" to screen
void Cat::Meow()
{
  cout << "Meow.\n";
}

// create a cat, set its age, have it
// meow, tell us its age, then meow again.
int main()
{
  Cat Frisky(5);
  Frisky.Meow();
  cout << "Frisky is a cat who is " ;
  cout << Frisky.GetAge() << " years old.\n";
  Frisky.Meow();
  Frisky.SetAge(7);
  cout << "Now Frisky is " ;
  cout << Frisky.GetAge() << " years old.\n";
  return 0;
}

Output: Meow.
Frisky is a cat who is 5 years old.
Meow.
Now Frisky is 7 years old.
