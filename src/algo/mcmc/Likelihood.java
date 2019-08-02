package algo.mcmc;



import java.io.IOException;

/**
 * Interface for Likelyhood
 * Distribution function for the logarithm of the likelihood 
 * @author fraison
 *
 */
public interface Likelihood {
	/**
	 * logarithm  of the likelihood value
	 * @param x
	 * @param psim
	 * @return
	 * @throws GaiaException
	 * @throws IOException
	 */
	   public double getLogLikelyhood(Parameters x, ParameterSimulation psim) throws IOException;
}
