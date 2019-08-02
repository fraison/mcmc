package algo.mcmc.test;

import java.io.IOException;

import algo.mcmc.Likelihood;
import algo.mcmc.ParameterSimulation;
import algo.mcmc.Parameters;


/**
 * Distribution function for the likelihood 
 * for the test here, it is actually a log likelihood
 * @author fraison
 *
 */

public class Likelihood4test implements Likelihood {
	double[] mu;
	double[] sigma;
	
	//@Override
	public double getLogLikelyhood(Parameters x, ParameterSimulation psim)
			throws IOException {
		
		return logMultiDimGaussian(x);
	}
	
	/**
	 * logarithm of a gaussian distribution
	 * @param x
	 * @return
	 */
	public  double logGaussian(Parameters x){ 
		   
		double cst = Math.sqrt(2. * Math.PI) ;
    	double a=0.;
        
    	a = 1./(cst*this.sigma[0])*Math.exp(-((x.getParameter(0)-this.mu[0])*(x.getParameter(0)-this.mu[0]))/(2.*this.sigma[0]*this.sigma[0])); 
    	return Math.log(a);
    
	}
 

	/**
	 * logarithm of a multi gaussian distribution
	 * @param x
	 * @return
	 */
	public double logMultiDimGaussian(Parameters x){ 
       
        double cst = -Math.log(2. * Math.PI) / 2.;
        
    	double a=0.;

        
    	for (int i =0;i<x.getParamSize();i++) a = a -((x.getParameter(i)-this.mu[i])*(x.getParameter(i)-this.mu[i]))/(2.*this.sigma[i]*this.sigma[i]) - Math.log(this.sigma[i]) + cst; 
    	return a;
    }
	

 
	/**
	 * setter of the parameters of the distribution for tests
	 * @param mu
	 * @param sigma
	 */
	public void setGaussian(double[] mu, double[] sigma){
		this.mu = mu.clone();
		this.sigma = sigma.clone();
	}
}
