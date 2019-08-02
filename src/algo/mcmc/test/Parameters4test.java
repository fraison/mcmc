package algo.mcmc.test;

import algo.mcmc.Parameters;

public class Parameters4test extends Parameters{

	public Parameters4test(double[] paramArray, double[] stepArray,
			boolean[] varying) {
		super(paramArray, stepArray, varying);
		
	}
	public  Parameters4test(int nparameters){
		super(nparameters);
	}
	
	@Override
	public Parameters deepCopy() {
		Parameters pm = new Parameters4test(this.param, this.step, this.vary);
		return pm;
		
	}

}
