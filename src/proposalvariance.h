

class ProposalVariance {

  public:
    // Constructor
    virtual ProposalVariance(double initial_pv,       // proposal variance
                             int adjust_iter, // adjust pv on multiples of adjust_iter
                             int max_iters);  // maximum iteration to adjust pv
    virtual ~ProposalVariance(); // Destructor
    void addaccept();  // Add to acceptance count
    void addreject();  // Add to iters but not accept count
    double getratio(); // 
    void resetratio(); 
    double getpv(); 
    void adjustpv();

  private:
    int accept_ct; // acceptance count
    int iter_ct;   // iteration count
    double pv;     // proposal variance
    double ratio;  // acceptance ratio (accept_ct/iter_ct)
    int adjust_iter; // iteration to adjust on
    int max_iter;

};

// ProposalVariance constructor
ProposalVariance::ProposalVariance(double initial_pv,
                                   int adjust_at_iter,
                                   int max_iters = 25000) {

  pv = initial_pv;
  adjust_at_iter = adjust_at_iter;
  max_iters = max_iters; 

}

// destructor
ProposalVariance::~ProposalVariance() { }


