package algo.mcmc.test;

import java.io.IOException;
import algo.mcmc.Likelihood;
import algo.mcmc.ParameterSimulation;
import algo.mcmc.Parameters;

public class Likelihood4test2  implements Likelihood {
	
	
	public static double logMultiPeakedGaussian(Parameters x){ 
		   
		double cst = Math.sqrt(2. * Math.PI) ;
    	double a=0.;
    	double[] mu = new double[]{100.,130.,600.};
    	double[] sigma = new double[]{30.,2.,7.};
    	double[] f = new double[]{1.,20.,1.};
        
    	for (int i =0;i<mu.length;i++) a = a +1./(cst*f[i]*sigma[i])*Math.exp(-((x.getParameter(0)-mu[i])*(x.getParameter(0)-mu[i]))/(2.*sigma[i]*sigma[i])); 
    	return Math.log(a);
    
    
 
	}

	public double getLogLikelyhood(Parameters x, ParameterSimulation psim)
			throws  IOException {
		
		return logMultiPeakedGaussian(x);
	}
}
