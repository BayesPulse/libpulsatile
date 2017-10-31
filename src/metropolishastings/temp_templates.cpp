#include <iostream>
#include <string>


template <typename writer>
class logger : public writer {
  public:
    void warn(const std::string& msg){ this->write("WARN", msg); }
    void error(const std::string& msg){ this->write("ERROR", msg); }
    void info(const std::string& msg){ this->write("INFO", msg); }
};


class file_writer {
  /*...*/
  public:
    void init( const std::string& file ) { /*...*/ }

  protected:
    void write( const std::string& kind, const std::string& msg ) {
      /*...*/
    }
};


class network_writer {
  /*...*/
  public:
    void init( const std::string& address, unsigned short port ) { /*...*/ }

  protected:
    void write( const std::string& kind, const std::string& msg ) {
      /*...*/
    }
};

typedef logger<file_writer> file_logger;
typedef logger<network_writer> network_logger;


int main() {

  network_logger nl;
  nl.init( "10.0.0.32", 8888 );
  nl.warn("this is another warning message");

}
