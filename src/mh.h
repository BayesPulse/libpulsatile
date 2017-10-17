// Abstract class for Metropolis Hastings samplers
// Demonstrates declaration of a constructors and destructor for the Cat class

//#include <iostream.h> // for cout
#include <Rcpp.h>
using namespace Rcpp;


class ModifiedMetropolisHastings  // begin declaration of the class
{

  public:                    // begin public section
    //Cat(int initialAge);   // constructor
    //~Cat();                // destructor
    virtual ModifiedMetropolisHastings(double starting_value,  // constructor
                                       double proposal_stddev, 
                                       Posterior posterior_function, 
                                       TruncLogic truncation,
                                       MCMC_state current_state); // all changing parameters?
    virtual ~ModifiedMetropolisHastings();           // destructor
    double sample(); // sample from posterior (run 1 iteration)
                     // Inputs/Outputs for sample()?


  private:                   // begin private section
    virtual posterior_function(); // logrho_calculation
    virtual parameter_support();         // truncation logic
    Proposalvariance pv;        // this should be an object of this class

};


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
