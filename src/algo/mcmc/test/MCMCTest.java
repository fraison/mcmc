package algo.mcmc.test;



import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.ResourceBundle;
import java.util.logging.Logger;
import javax.naming.ConfigurationException;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Test;
import algo.mcmc.ChainXconfig;
import algo.mcmc.Likelihood;
import algo.mcmc.MCMC;
import algo.mcmc.ParameterSimulation;
import algo.mcmc.Parameters;
import algo.mcmc.Prior;
import algo.mcmc.Proposal;
import fileUtils.Writer;
//import junit.framework;
/**
 * test for MCMC methods
 * @author fraison
 *
 */

public class MCMCTest  {
	private static final Logger LOGGER = Logger.getLogger( "MCMCTest");
	
	
	public void setUp() throws ConfigurationException {
		//ResourceBundle resource = ResourceBundle.getBundle("conf/mcmcTest");
		//PropertyLoader.load("mcmcTest.properties"); 
	}
	
	/**
	 * check that this QProp implementation for test (Proposal4test)  is providing a gaussian distributed samples
	 */
	@Test
	public void testqprop(){
		
		
		//parameter
		double tolerance = 0.1;
		double refmean = 12.;
		double refstd = 8.;
		double [] paramArray = new double[]{refmean};//starting value % real value = 10
		double [] stepArray = new double[]{refstd};// which is sigma*10
		boolean [] varying = new boolean[]{true};
		Parameters xt = new Parameters4test(paramArray, stepArray, varying);
		
		int niter = 10000;
		DescriptiveStatistics stats = new DescriptiveStatistics();
		
		for (int i=0;i<niter;i++){
			Proposal prop = new Proposal4test();
			Parameters y =  prop.getQProp(xt);
			stats.addValue(y.getParameter(0));
		}

		double mean = stats.getMean();		
		double std = stats.getStandardDeviation();		
		double skew = stats.getSkewness();		
		double kurt = stats.getKurtosis();
		
		LOGGER.info("sample mean="+mean+ "std="+std+" skew="+skew+ " kurt="+kurt);	
		
		// check that the collection of sample corresponds to the initial distribution 
		Assert.assertEquals(refmean , mean, tolerance);
		Assert.assertEquals(refstd , std, tolerance);
		Assert.assertEquals(0. , skew, tolerance);
		Assert.assertEquals(0. , kurt, tolerance);
		
		/*try {
			samplePlot(samples);
			sampleHistPlot(samples);
		} catch (GaiaException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}*/

		LOGGER.info("testqprop done");
		
	}
	
	/**
	 * sample the likelihood PDF (1D Gaussian here) at regularly separated points 
	 * and check that we get the expected max and CDF (at 3 sigma) of a Gaussian
	 * when prior is constant 
	 * @throws GaiaException
	 * @throws IOException
	 */
	@Test
	public void testlkd() throws Exception, IOException{
		double tolerance = 0.01;
		
		Proposal prop = new Proposal4test();
		Prior prior = new Pior4test();
		Likelihood lkd = new Likelihood4test();
		
		ParameterSimulation psim = new ParameterSimulation4test( prop, prior, lkd);
		
		//parameters 		
		double[] refmean = new double[]{ 20.};
		double[] refstd = new double[]{1.};

		//set likelihood for test
		((Likelihood4test) lkd).setGaussian(refmean, refstd);
		
		double [] paramArray = refmean.clone();//starting value % real value = 10
		double [] stepArray = refstd.clone();// which is sigma*10
		boolean [] varying = new boolean[]{true};
		
		Parameters xt = new Parameters4test(paramArray, stepArray, varying);
		
		int n= 100;//number of samples
		double[] lkdx = new double[n];
		double param;
		double cdf = 0.;//cumulative density function
		
		double span = 2* (3.*refstd[0]);//3 sigma on each side of the mean
		double min = refmean[0] - span/2.;
		for (int i=0;i<n;i++){
			param = min + i*span/n;
			xt.setParameter(0, param);
			lkdx[i] = Math.exp(lkd.getLogLikelyhood(xt, psim));//exponentiation to convert the log-likelihood back to likelihood 
			cdf = cdf+ lkdx[i]*span/n;
			
		}
		
		LOGGER.info("max="+lkdx[50]);
		LOGGER.info("cdf="+cdf);
		
		double max = 1./(Math.sqrt(2*Math.PI)*refstd[0]);//maximum of the distribution
		Assert.assertEquals(max , lkdx[50], tolerance);//check the maximum of the likelihood corresponds to the calculation for the mean
		Assert.assertEquals(1. , cdf, tolerance);//check that the cdf is close to 1 at 3 sigma
		
		//samplePlot(lkdx);

		LOGGER.info("testlkd done");
	}

	
	/**
	 * test the implementation of the logit function
	 */
	@Test
	public void testlogit(){
		
		double l1 = MCMC.logit(0.);
		LOGGER.info( "logit(0.) = "+ l1);
		
		double l2 = MCMC.logit(1.);
		LOGGER.info("logit(1.) = "+ l2);
		
		double l3 = MCMC.logit(0.5);
		LOGGER.info("logit(0.5) = "+ l3);
		
		
		Assert.assertEquals(Double.NEGATIVE_INFINITY,l1, 0D);
		Assert.assertEquals(Double.POSITIVE_INFINITY,l2, 0D);
		Assert.assertEquals(0.,l3, 0D);
		
	}

	/**
	 * Test the core function of MCMC: 1 basic MCMC chain sample return
	 * No loop here: we select a step value equal to 0. 
	 * Consequently the proposed new value will be the same as the old one. 
	 * Then the log likelihood of of old value and new will be identical. The acceptance rate 
	 * will be necessary 1 if the prior is constant
	 * @throws GaiaException
	 * @throws IOException
	 */
	@Test
	public void testptmcmc() throws IOException{
		double tolerance = 0.01;

		//parameters to recover			
		double[] refmean = new double[]{ 20.};
		double[] refstd = new double[]{2.};
		
		//starting values
		double [] paramArray = new double[]{20.};//starting value % real value = 10
		double [] stepArray = new double[]{0.};// which is sigma*10
		boolean [] varying = new boolean[]{true};
		Parameters xt = new Parameters4test(paramArray, stepArray, varying);
		
		//chain exponent (related to tempering): 1. for normal use
		double betai = 1.;
		
		
		Proposal prop = new Proposal4test();
		Prior prior = new Pior4test();
		Likelihood lkd = new Likelihood4test();
		//set likelihood for test
		((Likelihood4test) lkd).setGaussian(refmean, refstd);
		
		ParameterSimulation psim = new ParameterSimulation4test( prop, prior,  lkd);
		
		double  acc ;
		double logLkdy;
		
		//BASIC MCMCS	
		
		Object[] res;
		try {
			res = MCMC.ptmcmc(xt, betai, psim);		
			Parameters y = (Parameters4test)res[0];///in loop REUSE XT (deep copy)
			acc = (Double) res[1];
			logLkdy= (Double) res[2];
			
			//acceptance rate is 1 
			Assert.assertEquals(1. , acc, 0.);
			
			//new parameter has the same value as the initial one (20)
			Assert.assertEquals(xt.getParameter(0) , y.getParameter(0), 0.);
			
			//check that the likelihood corresponds to the max of the gaussian for step =0
			Assert.assertEquals(Math.log(1./(Math.sqrt(2*Math.PI)*refstd[0]) ), logLkdy, tolerance);
			LOGGER.info("testptmcmc done");
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		
	}

	
	
	/**
	 * Test 1 basic MCMC chain: we sample a 1D Gaussian (mean=20, std=2)
	 * and start from 15 with a step =4. We do 1000 iterations (get rid of the  50 first values 
	 * of the burning period).
	 * In the proposal function, we generate new values with half the probability as step = 2 sigma
	 * and step is used as a new sigma for generation. This is why the acceptance rate is expected to be 1/2.
	 * @throws GaiaException
	 * @throws IOException
	 */
    @Test
	public  void testsimpleChainX() throws Exception, IOException{
		double tolerance = 0.5;
		double tolearnceRate = 0.05;
		//parameters to recover			
		double[] refmean = new double[]{ 20.};
		double[] refstd = new double[]{2.};
		
		//starting values
		double [] paramArray = new double[]{15.};//starting value % real value = 10
		double [] stepArray = new double[]{4.};// which is sigma*10
		boolean [] varying = new boolean[]{true};
		Parameters xt = new Parameters4test(paramArray, stepArray, varying);
		
		//chain exponent (related to tempering): 1. for normal use
		double betai = 1.;
		
		
		Proposal prop = new Proposal4test();
		Prior prior = new Pior4test();
		Likelihood lkd = new Likelihood4test();
		//set likelihood for test
		((Likelihood4test) lkd).setGaussian(refmean, refstd);
		
		ParameterSimulation psim = new ParameterSimulation4test( prop, prior,  lkd);
		int niter = 1000;
		double[] samples;
		double[]  acc ;
		
		//BASIC MCMCS		

    	Object [] res = MCMC.simpleChainX(niter, xt, betai, psim);
    	Parameters[] yp =  MCMC.getChainXParameters(res);
    	samples =  Parameters.extractor(yp,0);    	
    	acc = MCMC.getChainXAcceptance(res);
 	

    	
    	double accmean = StatUtils.mean(acc.clone());
		LOGGER.info("acceptance rate="+accmean);	
	
		int endBurn = 100; //burning period: to be taken out
	
		//GVectorNd data = new GVectorNd(Array1DUtilities.vecCopy(endBurn,samples));	
				
		DescriptiveStatistics stats = new DescriptiveStatistics();
		
		for (int i=endBurn;i<niter;i++){
			stats.addValue(samples[i]);
		}
        double median = stats.getPercentile(50.); //median
		double mean = stats.getMean();		
		double std = stats.getStandardDeviation();		
		double skew = stats.getSkewness();		
		double kurt = stats.getKurtosis();

		
		LOGGER.info("sample median="+median+" sample mean="+mean+ " std="+std+" skew="+skew+ " kurt="+kurt);
		
		//plots.... 
		//samplePlot(samples);
		//sampleHistPlot(samples);
		
		//check that the mean of the sample corresponds to the reference distribution
		Assert.assertEquals(refmean[0] , mean, tolerance);
		Assert.assertEquals(refstd[0] , std, tolerance);
		Assert.assertEquals(0. , skew, tolerance);
		Assert.assertEquals(0. , kurt, 10*tolerance);
		
		double refrate = 0.5;
		Assert.assertEquals(refrate, accmean, tolearnceRate);

		
		LOGGER.info("testsimpleChainX done");			
		
	}
    
    
	/**
	 * Test calculation of the optimal step
	 * We generate a set of steps according to a reference slope
	 * and check we can recover the optimal step parameter from the method 
	 * doing the calculation (no chain is running there)
	 * @throws GaiaException
	 * @throws IOException
	 */
	@Test
	public  void testStepCalculator2() throws  IOException{
		
		double tolerance = 1e-7;
		
		//we want to find the step to get an acceptance rate of 0.3 (30%)
		double  targetRate = 0.3;
		boolean diagnostic = true;
		
		//we start from a reference slope:
		double refSlope = -1.01;
		double refOffset = 1.2 ; 
		
		double startStep  = 4.;
		
		int trials = 11; //like in code
     		
		List<Double> mrate= new ArrayList<Double>() ;
		List<Double> mstep = new ArrayList<Double>() ;
		double lstep;
		double lcalcRate;
		
		//we generate the list of steps and acceptance rates
		for (int i = 0; i < trials; i++){
			mstep.add( startStep*10./Math.pow(2., i));  
			
			lstep = Math.log(mstep.get(i));
			lcalcRate= ((refSlope*lstep+ refOffset) );
			//this is the main formula connecting rate and step
			mrate.add( Math.exp(lcalcRate) / ( Math.exp(lcalcRate)+1.) );
		}
		
		//calculate the step value that the tested method should recover
		double refStep =  Math.exp((MCMC.logit(targetRate)-refOffset)/refSlope  ); 
			
			
		//double step 
		
		double[] sol;
		try {
			sol = MCMC.stepCalculator2(targetRate, mrate.toArray(), mstep.toArray(),diagnostic,0);
			double optistep = sol[0];
			double slope = sol[1];
			double offset = sol[2];
			
			LOGGER.info("optistep="+optistep+ " slope="+slope+ " offset+"+offset);		
			
			Assert.assertEquals(refSlope,slope, tolerance);
			Assert.assertEquals(refOffset,offset, tolerance);
			Assert.assertEquals(refStep,optistep, tolerance);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}//{optiStep,slope,offset}
		

	}
	
	/**
	 * Test optimal step determination
	 * 
	 * @throws GaiaException
	 * @throws IOException
	 */
	@Test
	public  void test1StepSelectorX2() throws  IOException{
		
		double tolearnceRate = 0.05;
		//parameters to recover			
		double[] refmean = new double[]{ 20.};
		double[] refstd = new double[]{2.};
		
		//starting values
		double [] paramArray = new double[]{15.};//starting value % real value = 10
		double [] stepArray = new double[]{4.};// lead to acceptance rate ~0.5 (50%) with value = 4.

		
		//TODO crazy things with STEpARRAy
		//3 test: 1 normal, 2 huge step, 3 tiny step and monitor step reduction
		
		
		boolean [] varying = new boolean[]{true};
		Parameters xt = new Parameters4test(paramArray, stepArray, varying);
		
		//chain exponent (related to tempering): 1. for normal use
		double betai = 1.;
		
		
		Proposal prop = new Proposal4test();
		Prior prior = new Pior4test();
		Likelihood lkd = new Likelihood4test();
		//set likelihood for test
		((Likelihood4test) lkd).setGaussian(refmean, refstd);
		
		ParameterSimulation psim = new ParameterSimulation4test( prop, prior,  lkd);
		int niter = 1000;
		double[]  acc ;
		
		
		//it will try to find the step to get an acceptance rate of approx. 0.3 (30%)
		double  targetRate = 0.3;
		

		boolean diagnostic = true;
		int minorCycles =100;//100 or 50: number of iteration for each individual chain
	    double gamma = 10; //very important here
    	int trials = 11;//11 number of chains to run : 1 per trial step
    	int ncycles =1;
    	
		
		String s2 = "./data/test/"; //should be under mcmc in order to work
		fileUtils.Writer w2 = MCMC.log2(1, trials, s2);


		double optistep;
		try {
			optistep = MCMC.stepSelectorX2(diagnostic,betai,targetRate, xt, minorCycles, gamma, w2, psim,trials, ncycles, 0,0)[0];
			LOGGER.info("optistep="+optistep);
			
			
			//Check that the calculated step leads to the target rate with a chain
			xt.setSteps(0, optistep);
			
			Object [] res = MCMC.simpleChainX(niter, xt, betai,psim);
			
			//Parameters[] yp =  MCMC.getChainXParameters(res);
			//double[] samples = Parameters.extractor(yp,0);
			
			acc = MCMC.getChainXAcceptance(res);		
			
		
			//results ...
	    	double accmean = StatUtils.mean(acc.clone());
			LOGGER.info("acceptance rate="+accmean);	
			

			Assert.assertEquals(targetRate, accmean, tolearnceRate);//functionality to test
		
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}


		
		//cleanup
		//String filename = w2.getFileName();
		w2.deletefile();


		LOGGER.info("test1StepSelectorX2 done");
			
		
	}

	
	
	/**
	 * Test optimal step selection when step much too small (rate very high)
	 * 
	 * @throws GaiaException
	 * @throws IOException
	 */
	@Test
	public  void test2StepSelectorX2() throws IOException{
		

		//parameters to recover			
		double[] refmean = new double[]{ 20.};
		double[] refstd = new double[]{2.};
		
		//starting values
		double [] paramArray = new double[]{15.};//starting value % real value = 10
		double [] stepArray = new double[]{0.04};// 4. leads to acceptance rate ~0.5 (50%) -> 0.04 way too small
		
		
		boolean [] varying = new boolean[]{true};
		Parameters xt = new Parameters4test(paramArray, stepArray, varying);
		
		//chain exponent (related to tempering): 1. for normal use
		double betai = 1.;
		
		
		Proposal prop = new Proposal4test();
		Prior prior = new Pior4test();
		Likelihood lkd = new Likelihood4test();
		//set likelihood for test
		((Likelihood4test) lkd).setGaussian(refmean, refstd);
		
		ParameterSimulation psim = new ParameterSimulation4test( prop, prior,  lkd);
		int niter = 1000;
		///double[] samples;
		double[]  acc ;
		
		
		//find the step to get an acceptance rate of 0.3 (30%)
		double  targetRate = 0.3;


		
		boolean diagnostic = true;
		int minorCycles =100;//100 50 number of iteration for each individual chain
	    double gamma = 10;
    	int trials = 11;//11 number of chains to run : 1 per trial step
    	int ncycles = 1;
    	
		String s2 = "data/test/";			
		Writer w2 = MCMC.log2(1, trials, s2);
		
	    //initial chain
		Object[] res;
		try {
			res = MCMC.simpleChainX(niter, xt, betai, psim);
			
			
			//Parameters[] yp =  MCMC.getChainXParameters(res);
			//samples =  Parameters.extractor(yp,0);
			
			acc = MCMC.getChainXAcceptance(res);		
			
		
			//results ...
	    	double accmeanInit = StatUtils.mean(acc.clone());////take clone as statistics modifies input
			LOGGER.info("initial acceptance rate="+accmeanInit);	
		    
		    
			//calculate new  step
			double optistep = MCMC.stepSelectorX2(diagnostic,betai,targetRate, xt, minorCycles, gamma, w2, psim, trials, ncycles, 0,0)[0];
			LOGGER.info("optistep="+optistep);
			
			
			//Check that the calculated step leads to the target rate 
			xt.setSteps(0, optistep);
			
			//chain
			res = MCMC.simpleChainX(niter, xt, betai, psim);
			
			//yp =  MCMC.getChainXParameters(res);
			//samples =  Parameters.extractor(yp,0);
			
			acc = MCMC.getChainXAcceptance(res);		
			
		
			//results ...
	    	double accmean = StatUtils.mean(acc.clone());////take clone as statistics modifies input
			LOGGER.info("final acceptance rate="+accmean);	
					
		
			Assert.assertTrue("reduce acceptance rate", (accmean<accmeanInit));

						
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		//cleanup
		String filename = w2.getFileName();
		w2.deletefile();
		
		LOGGER.info("test2StepSelectorX2 done");
			
		
	}

	/**
	 * Test 1 MCMC chain with automatic optimal step calculation with 2 parameters
	 * @throws GaiaException
	 * @throws IOException
	 */
	@Test
	public void testChainX() throws IOException{
		double tolerance = 0.1;
		double tolearnceRate = 0.1;
		//parameters to recover			
		double[] refmean = new double[]{ 20.,100};
		double[] refstd = new double[]{2.,10.};
		
		//starting values
		double [] paramArray = new double[]{15.,110.};//starting value % real value = 10
		double [] stepArray = new double[]{4.,8.};// which is sigma*10
		boolean [] varying = new boolean[]{true, true};
		Parameters xt = new Parameters4test(paramArray, stepArray, varying);
		
		//chain exponent (related to tempering): 1. for normal use
		double betai = 1.;
		
		
		Proposal prop = new Proposal4test();
		Prior prior = new Pior4test();
		Likelihood lkd = new Likelihood4test();//multi dimensional gaussian
		//set likelihood for test
		((Likelihood4test) lkd).setGaussian(refmean, refstd);
		
		ParameterSimulation psim = new ParameterSimulation4test( prop, prior,  lkd);
		
		double[] samples0;
		double[] samples1;
		double[]  acc ;

		int trials = 11;
		
		//MCMCS	with automatic optimal step tuning
		double  targetRate = 0.30;
		
		String s2 = "data/test/";	
		Writer w2 = MCMC.log2(1, trials, s2);
		
		
		//get configuration of the run from properties
		ChainXconfig cfg;
		try {
			//cfg = new ChainXconfig("conf/mcmcTest.properties");
			cfg = new ChainXconfig("conf/mcmcTest");
			//cfg = new ChainXconfig("/home/fraison/rzgdatashare/pro/softEng/workspaceMCMC/mcmc/src/conf/mcmcTest");
			cfg.targetRate = 0.30;		
			cfg.iter =10000;
			cfg.diagnostic = true;
			
			
			//run the chain
			Object [] res = MCMC.chainX(xt, betai, psim, w2 ,cfg);
			
			Parameters[] yp =  MCMC.getChainXParameters(res);
			samples0 =  Parameters.extractor(yp,0);//param 0
			samples1 =  Parameters.extractor(yp,1);// param 1
			double[] loglkd = (double[])res[2];
			
			acc = MCMC.getChainXAcceptance(res);		
		
			//results ...
	    	double accmean = StatUtils.mean(acc.clone());////take clone as statistics modifies input
			LOGGER.info("final acceptance rate="+accmean);	
			LOGGER.info("acceptance rate="+accmean);	
			
		
			//VALIDATION
		
			int endBurn = 100; //burning period: to be taken out
				
			DescriptiveStatistics stats = new DescriptiveStatistics();
			
			for (int i=endBurn;i<samples0.length;i++){
				stats.addValue(samples0[i]);
			}
	        double median0 = stats.getPercentile(50.); //median
			double mean0 = stats.getMean();		
			double std0 = stats.getStandardDeviation();		
			double skew0 = stats.getSkewness();		
			double kurt0 = stats.getKurtosis();
			
			LOGGER.info("sample median="+median0+" sample mean="+mean0+ " std="+std0+" skew="+skew0+ " kurt="+kurt0);	
			
			
			DescriptiveStatistics stats1 = new DescriptiveStatistics();
			
			for (int i=endBurn;i<samples1.length;i++){
				stats1.addValue(samples1[i]);
			}
	        double median1 = stats1.getPercentile(50.); //median
			double mean1 = stats1.getMean();		
			double std1 = stats1.getStandardDeviation();		
			double skew1 = stats1.getSkewness();		
			double kurt1 = stats1.getKurtosis();
			LOGGER.info("sample median="+median1+" sample mean="+mean1+ " std="+std1+" skew="+skew1+ " kurt="+kurt1);			//plots.... 
			
			/*samplePlot(samples0);
			sampleHistPlot(samples0);
			samplePlot(samples1);
			sampleHistPlot(samples1);*/
			
			Assert.assertEquals(targetRate, accmean, tolearnceRate);
			Assert.assertEquals(refmean[0] , mean0, tolerance);
			Assert.assertEquals(refstd[0] , std0, tolerance);
			Assert.assertEquals(0. , skew0, tolerance);
			Assert.assertEquals(0. , kurt0, 10*tolerance);
			
			
			LOGGER.info("sample MAP param 0 ="+MCMC.getSamplesMap(0, 0., 1000., samples0, loglkd));
			LOGGER.info("sample MAP param 1 ="+MCMC.getSamplesMap(0, 0., 1000., samples1, loglkd));

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
		
		//cleanup
		String filename = w2.getFileName();
		w2.deletefile();
		
		LOGGER.info("testChainX done");
			
		
	}
	
	/**
	 * test chain swapping test
	 * TODO Check if the logic of the test is fine here if test doesn't execute
	 * @throws GaiaException
	 * @throws IOException
	 */
	@Test
	public void testptAccept() throws IOException{		
		//parameters to recover			
		double[] refmean = new double[]{ 20.};
		double[] refstd = new double[]{2.};
		
		Proposal prop = new Proposal4test();
		Prior prior = new Pior4test();
		Likelihood lkd = new Likelihood4test();//multi dimensional gaussian
		
		//set likelihood for test
		((Likelihood4test) lkd).setGaussian(refmean, refstd);
		ParameterSimulation psim = new ParameterSimulation4test( prop, prior,  lkd);
		
		//TEST1
		//starting values
		double [] paramArray = new double[]{15.};//starting value % real value = 10
		double [] stepArray = new double[]{0.04};// 4. leads to acceptance rate ~0.5 (50%) -> 0.04 way too small			
		boolean [] varying = new boolean[]{true};
		Parameters xti = new Parameters4test(paramArray, stepArray, varying);
		
		
		double [] paramArray1 = new double[]{20.};//starting value % real value = 10
		double [] stepArray1 = new double[]{0.05};// 4. leads to acceptance rate ~0.5 (50%) -> 0.04 way too small			
		boolean [] varying1 = new boolean[]{true};
		Parameters xti1 = new Parameters4test(paramArray1, stepArray1, varying1);
		
		
		double betai = 0.5;
		double betai1 = 0.5;
		boolean try1;
		
		boolean test1 = true;
		for (int i =0;i<1100;i++){
			try {
				try1 = MCMC.ptAccept(xti, betai, xti1, betai1, psim);
				if (!try1 ) test1 = false;
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		}
		
		//TEST2
		//starting values
		paramArray = new double[]{15.};//starting value % real value = 10
		stepArray = new double[]{0.04};// 4. leads to acceptance rate ~0.5 (50%) -> 0.04 way too small			
		varying = new boolean[]{true};
		xti = new Parameters4test(paramArray, stepArray, varying);
		
		
		paramArray1 = new double[]{15.};//starting value % real value = 10
		stepArray1 = new double[]{0.04};// 4. leads to acceptance rate ~0.5 (50%) -> 0.04 way too small			
		varying1 = new boolean[]{true};
		xti1 = new Parameters4test(paramArray1, stepArray1, varying1);
		
		betai = 0.5;
		betai1 = 0.7;
		
		boolean test2 = true;
		for (int i =0;i<1100;i++){
			try {
				try1 = MCMC.ptAccept(xti, betai, xti1, betai1, psim);
				if (!try1 ) test2 = false;
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		}

		Assert.assertEquals(true, test1);
		Assert.assertEquals(true, test2);
	}
	
	
	/**
	 * test parallel tempering on a multimodal 1D distribution
	 * @throws GaiaException
	 * @throws IOException
	 */
	@Test
	public void testptChainX() throws Exception, IOException{
		
		ResourceBundle resource = ResourceBundle.getBundle("conf/mcmcTest");
		int niter = 10000;
		//get configuration of the run from properties
		ChainXconfig cfg = new ChainXconfig("conf/mcmcTest");	
		double [] betai = cfg.betai;
		
		Proposal prop = new Proposal4test();
		Prior prior = new Pior4test();
		Likelihood lkd = (Likelihood) new Likelihood4test2();//multipeaked gaussian		
		ParameterSimulation psim = new ParameterSimulation4test( prop, prior,  lkd);
		
		//starting values (only 1 parameter used here)
		double [] paramArray = new double[]{170.};
		double [] stepArray = new double[]{2.};
		boolean [] varying = new boolean[]{true};
		Parameters xin = new Parameters4test(paramArray, stepArray, varying);
		
		final MCMC  mymcmc = new MCMC();
		mymcmc.ptChainX( xin,psim, cfg);
		
    	int psize = paramArray.length;    	
    	int lastChain = betai.length-1; 
    	
		for(int n =0;n<psize;n++){
			
			
			double[] y =new double[niter] ;

			double[] loglkd =  mymcmc.getChainLogLikelihoods()[lastChain];
			double[] samples = mymcmc.getSamples(lastChain, 0);
			 
			double map0 = MCMC.getSamplesMap(lastChain, 0., 110., samples, loglkd);
			double map1 = MCMC.getSamplesMap(lastChain, 110., 200., samples, loglkd);
			double map2 = MCMC.getSamplesMap(lastChain, 400, 800., samples, loglkd);
			 
			LOGGER.info("sample MAP0 param 0 ="+map0);
			LOGGER.info("sample MAP1 param 0 ="+map1);
			LOGGER.info("sample MAP2 param 0 ="+map2);
			
			double[] mu = new double[]{100.,130.,600.};
			double tolerance =1;
			Assert.assertEquals(mu[0], map0, tolerance);
			Assert.assertEquals(mu[1], map1, tolerance);
			Assert.assertEquals(mu[2], map2, tolerance);
			


		}
		cfg.saveXconfig("");
		
		//cleanup
		//String filename = w2.getFileName();
		//FileUtilities.deletefile(filename);
		
		LOGGER.info("testptChainX done");
		
		
	}


	
	public static void main (String[] args) {
		 try {
			new MCMCTest().testqprop();
			//new MCMCTest().testsimpleChainX();
			//new MCMCTest().test2StepSelectorX2();
			//new MCMCTest().testChainX();
			//new MCMCTest().testptChainX();
			// LOGGER.info("test 1");
			// new MCMCTest().setUp();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
}
