
#include <armadillo>
//#include "catch.h"
#include "counter.h"
#include "proposalvariance.h"

int main() {

  std::cout << "Hello world\n";


  ////////////////////////////////////////
  // testing the weird issue counter class having wrong result of getratio()
  //  (the loop and counting here works right, so something is wrong in the
  //  class)
  ////////////////////////////////////////
  int i;
  int accept = 0;
  int iter = 0;
  double proportion = 0.0;

  for (int i = 0; i < 100; ++i) {
    if (i % 4 == 0) accept++;
    iter++;

    proportion = double(accept) / double(iter);
    std::cout << "Num accepted is: " << accept << " Num rejected is: " 
      << iter << " With a ratio of: " << proportion << "\n";
  }
  ////////////////////////////////////////

  return 0;

}


