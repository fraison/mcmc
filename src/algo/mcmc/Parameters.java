package algo.mcmc;



/**
 * Class with the main methods for MCMC to deal with parameters. For MCMC, each input parameter should have an associated step value 
 * and a boolean telling if the parameter is varying or not.
 * A constructor with the size of the parameter array is accessible in case the parameter needs to be used outside of the MCMC algorithm. But Dangerous.
 * The class has to be extended to create the mapping between the list of parameters (double[]) and their meaning for the specific problem which is solved.
 * @author fraison
 *
 */

public abstract class Parameters /*extends simpleParameter (no steps neither varying) ? */ {
	int nparam;
	protected double[] param;
	protected double[] step;
	protected boolean[] vary;
	
	/**
	 * constructor for MCMC
	 * @param paramArray
	 * @param stepArray
	 * @param varying
	 */
	public Parameters(double[] paramArray, double[] stepArray, boolean[] varying){
		this.param = paramArray.clone();
		if (stepArray!=null) this.step = stepArray.clone();
		if (varying!=null) this.vary = varying.clone();
		this.nparam = this.param.length;		
		

	}
	
	/**
	 * Constructor out of MCMC for simulations for example
	 * @param nparameters
	 */
	public  Parameters(int nparameters){
		this.nparam = nparameters;
		this.param = new double [this.nparam];
		this.step = new double [this.nparam];
		this.vary = new boolean [this.nparam];
	}
	
	
	/**
	 * Getter for parameters
	 * @return
	 */
	public double []getParameters(){
		return this.param;
	}

	/**
	 * 	Getter for a parameter
	 * @param i index in array
	 * @return parameter
	 */
	public double getParameter(int i){
		return this.param[i];
		
	}
	
	/**
	 * Setter for parameters
	 * @param paramArray array of parameters
	 */
	public void setParameters(double[] paramArray){
		this.param = paramArray.clone();
		this.nparam = this.param.length;
	}
	
	/**
	 * Setter for a parameter
	 * @param i index in array
	 * @param param
	 */
	public void setParameter(int i, double param){
		this.param[i] = param;
		
	}
	
	/**
	 * Getter for steps
	 * @return
	 */
	public double []getSteps(){
		return this.step;
	}
	
	/**
	 * Getter for a step
	 * @param i index in array
	 * @return
	 */
	public double getStep(int i){
		return this.step[i];
		
	}
	
	/**
	 * Setter for steps
	 * @param stepArray
	 */
	public void setSteps(double[] stepArray){
		this.step = stepArray.clone();
		
	}
	
	/**
	 * Setter for a step
	 * @param i index in array
	 * @param step
	 */
	public void setSteps(int i, double step){
		this.step[i] = step;
		
	}
	
	/**
	 * get Varying array
	 * @return array of varying parameter
	 */
	public boolean[] getVaryings(){
		return this.vary;
	}
	
	/**
	 * Getter for varying
	 * @param i index in array
	 * @return
	 */
	public boolean getVarying(int i){
		return this.vary[i];
	}
	
	/**
	 * Setter for varying
	 * @param v
	 */
	public void setVaryings(boolean[] v){
		this.vary = v.clone();
	}
	
	/**
	 * Setter for varying
	 * @param i index in array
	 * @param v
	 */
	public void setVaryings(int i, boolean v){
		this.vary[i] = v;
	}

	/**
	 * Size of the parameter array 
	 * @return number of parameters
	 */
	public int getParamSize(){
		return nparam;
	}
	
	/**
	 * Deep copy but cannot use constructor
	 * @return
	 */
	public abstract Parameters deepCopy();

	/**
	 * Useful?
	 * @param n
	 * @return
	 */
	public static int swapper(int n){

		if (n==0) return 1;
		if (n==1) return 0;
		else throw new IllegalArgumentException("problem with inital step ouput");
			
	}
	
	/**
	 * Useful?
	 * extract from a list of parameters, the values corresponding to a single parameter for each iteration
	 * @param yl array of vector of parameters
	 * @param p order of the requested parameter in the vector
	 * @return array of value of parameter number p
	 */
	public static double[] extractor(Parameters[] yl, int p){
		int iter  = yl.length;
		double[] paramp = new double[iter];
		for (int i=0;i<iter;i++)     		paramp[i] = yl[i].getParameter(p);    	
		return paramp;
	}



	

}
