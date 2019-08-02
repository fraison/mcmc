package algo.mcmc;


/**
 * 
 * top class for the proposal function (transition kernel), Prior and likelihood
 * @author fraison
 *
 */
public abstract class ParameterSimulation {
	
	Proposal lprop;
	Prior lprior;
	Likelihood llkd;
	
	public ParameterSimulation(Proposal prop, Prior prior, Likelihood lkd ) {
		this.lprop = prop;
		this.lprior = prior;
		this.llkd = lkd;
		
	}
	
	public Proposal getProposalClass(){
		return lprop;
	}
	public Prior getPriorClass(){
		return lprior;
	}
	public Likelihood getLikelihoodClass(){
		return llkd;
	}

}
