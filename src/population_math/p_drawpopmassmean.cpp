
class drawPopMassMean : public Gibbs
{

  private:
    bool parameter_support();
    double posterior_function(Patient patient);
    JointProposalVariance pv();

}


bool parameter_support(double proposal)
{

};

