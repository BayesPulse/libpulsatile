

#include <testthat.h>
#include "counter.h"

//
// Test the Counter class
//

int main() {
  printf("Counter can count\n");

  Counter cnt(0.25, 500, 25000);

  cnt.addreject();
  printf("Expecting 1", cnt.getiter());

  cnt.addaccept();
  expect_true( cnt.getaccept() == 1  );
  expect_true( cnt.getiter() == 2    );
  expect_true( cnt.getratio() == 0.5 );

  cnt.resetratio();
  expect_true( cnt.getaccept() == 0  );
  expect_true( cnt.getiter() == 0    );
  expect_true( cnt.getratio() != cnt.getratio() );

  }
}

