

int main() {



};


struct myContainer
{

  double a, b;
  int c;

  // Constructor
  myContainer(ia, ib, ic) {
    a = ia;
    b = ib;
    c = ic;
  };

}


class operatorClass
{

  operatorClass() {
    // How to on implementation of class, tell the function
    // operate_on_list_of_containers() which parameter to operate on?
  };

  double operate_on_list_of_containers(std::list<myContainer> mylist) {
    std::list<myContainer>::iterator iter      = mylist.begin();
    std::list<myContainer>::const_iterator end = mylist.end();

    double sum_of_squares = 0;
    while ( iter != end ) {
      sum_of_squares += (iter->a * iter->a);
    }

    return sum_of_squares;

  };

}
