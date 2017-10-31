#include <iostream>
#include <string>
#include <vector>

class logger {
  public:
    void warn(const std::string& msg){ write("WARN", msg); }
    void error(const std::string& msg){ write("ERROR", msg); }
    void info(const std::string& msg){ write("INFO", msg); }

    virtual ~logger() { }

  protected:
    virtual void write( const std::string& kind, const std::string& msg ) = 0;
};

class file_logger: public logger {
  /*...*/
  public:
    file_logger( const std::string& file ) { /*...*/ }
    ~file_logger() { /*...*/ }

  protected:
    virtual void write( const std::string& kind, const std::string& msg ) {
      /*...*/
    }
};

class network_logger: public logger {
  /*...*/
  public:
    network_logger( const std::string& address, unsigned short port ) { /*...*/ }
    ~network_logger() { /*...*/ }

  protected:
    virtual void write( const std::string& kind, const std::string& msg ) {
      /*...*/
    }
};



int main() {

  // as pointer to the base class
  logger *l = new file_logger( "foo.log" );
  l->warn("this is a warning message");

  // as automatic variable on the stack
  network_logger nl( "10.0.0.32", 8888 );
  nl.warn("this is another warning message");

}


