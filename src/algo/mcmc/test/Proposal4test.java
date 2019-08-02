package algo.mcmc.test;

import algo.mcmc.MCMC;
import algo.mcmc.Parameters;
import algo.mcmc.Proposal;

public class Proposal4test implements Proposal{
     /**
      * Proposal function (transition kernel) : Metropolis-Hasting (Gaussian generator)
      * It provides a new value for each parameter
      */
	//@Override
	public Parameters getQProp(Parameters x) {
		double inter;
    	Parameters xout = new Parameters4test(x.getParamSize());
    	for (int i = 0; i<x.getParamSize();i++){
    		inter = MCMC.gaussGen( x.getParameter(i), x.getStep(i));
    		xout.setParameter(i, inter);
    		xout.setSteps(i, x.getStep(i));
    		xout.setVaryings(i, x.getVarying(i));
    	}
		return xout;
	}

}
