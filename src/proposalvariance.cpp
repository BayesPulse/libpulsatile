
// Constructor
void ProposalVariance::ProposalVariance(double initial_pv,
                                        int adjust_at_iter,
                                        int max_iters = 25000); {

  pv             = initial_pv;
  adjust_at_iter = adjust_at_iter;
  max_iters      = max_iters;
  accept_ct      = 0;
  iter_ct        = 0;
  ratio          = 0;

}

// destructor -- need to define here?
ProposalVariance::~ProposalVariance() { }

void ProposalVariance::adjustpv()
{

}
