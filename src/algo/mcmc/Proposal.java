package algo.mcmc;



/**
 * interface for a transition kernel
 * @author fraison
 *
 */
public interface Proposal {
	/**
	 * this function should provide a new proposed array of parameter values
	 **/
	public  Parameters getQProp( Parameters x);
}
