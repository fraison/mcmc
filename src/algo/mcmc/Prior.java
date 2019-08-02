package algo.mcmc;




/**
 * Interface of Priors
 * @author fraison
 *
 */
public interface Prior {
	
	/**
	 * Return the product of all prior values 
	 * @param xt parameters
	 * @param psim additional simulation data
	 * @return
	 */
	public  double getPrior(Parameters xt, ParameterSimulation psim);
}
