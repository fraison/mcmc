package algo.mcmc;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;
import java.util.Random;
import java.util.logging.Logger;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import fileUtils.Writer;
/**
 * Monte Carlo Markov Chain solver
 * @author fraison
 * TODO cleanup
 */

public class MCMC {
	private static final Logger LOGGER = Logger.getLogger( "MCMC");
	
    static Random rand = new Random(0L);  
    
    Object[] out ;
    double[] lbetai;
    
    /**
     * prepare the logging of step tuning phase (stepSelectorX2)
     * @param nbetai number of exponents 
     * @param nparam number of parameters
     * @param dirName directory name
     * @return
     * @throws FileNotFoundException
     */
    public static Writer log1(int nbetai,int nparam, String dirName) throws FileNotFoundException{
    	String s1 = dirName+"log1_" + ( new SimpleDateFormat ("yyMMdd-HHmmss") ).format(Calendar.getInstance().getTime())+".dat";
    	Writer w1 = new Writer(s1);
    	w1.addValFormat("iter", "%8d", "%8s");
    	w1.addValFormat("swapRate", "%10.3f", "%9s");
    	for (int b =0;b<nbetai;b++){
    		w1.addValFormat(("rate["+b+"]"),"%10.3f","%10s");
    		for (int l=0;l<nparam;l++)	w1.addValFormat(("oldstep["+b+"]["+l+"]"),"%14.3e", "%14s");
    		for (int l=0;l<nparam;l++)	w1.addValFormat(("step["+b+"]["+l+"]"),"%14.3e", "%14s");	
    	}
    	return w1;
    }
    
    
    /**
     * prepare the logging of step tuning phase (stepSelectorX2)
     * @param nbetai number of exponents 
     * @param trialsOptStep number of trials to interpolate optimal step
     * @param dirName directory name
     * @return
     * @throws FileNotFoundException
     */
    public static Writer log2(int nbetai,int trialsOptStep, String dirName) throws FileNotFoundException{
    	String s2 = dirName+"log2_" + ( new SimpleDateFormat ("yyMMdd-HHmmss") ).format(Calendar.getInstance().getTime())+".dat";	
    	Writer w2 = new Writer(s2);
    	LOGGER.info("temporary data log file:"+s2);
    	w2.addValFormat("iter", "%8d", "%8s");
    	w2.addValFormat("chain", "%6d", "%6s");
    	w2.addValFormat("param", "%6d", "%6s");
    	w2.addValFormat("cycle", "%6d", "%6s");
    	w2.addValFormat("type", "%5s", "%5s");
    	for (int l=0;l<trialsOptStep;l++) w2.addValFormat(("step["+l+"]"),"%11.3e", "%11s");
    	for (int l=0;l<trialsOptStep;l++) w2.addValFormat(("rate["+l+"]"),"%11.3e", "%11s");
    	w2.addValFormat("opti", "%11.3e", "%11s");
    	w2.addValFormat("mode", "%6s", "%6s"); 	
    	w2.addValFormat("slope", "%11.3e", "%11s");
    	w2.addValFormat("offset", "%11.3e", "%11s"); 	  	

    	return w2;
    	
    	
    }
    
    //TODO ? useful
    
    public  void ptChainX( Parameters xin, ParameterSimulation psim, ChainXconfig cfg) throws Exception, IOException{
   	
    	//control parameters

    	
    	double targetRate = cfg.targetRate;
    	double swapTargetRate = cfg.swapTargetRate;
    	
    	boolean optStep = cfg.optStep;
    	
    	int n1 = cfg.majorCycles;
    	int minorCycles = cfg.minorCycles;		

	    double lambda = cfg.lambda;
	    double lambdaScale = cfg.lambdaScale;
	    double gamma = cfg.gamma;
	    
    	int trialsOptStep = cfg.trialsOptStep;
    	int ncyclesOptStep = cfg.ncyclesOptStep;	
    	
    	int iter = cfg.iter;    	
    	double betai[] = cfg.betai;
    	lbetai = betai.clone();
    	
    	String   dirName = cfg.dirName;
    	
    	
    	Parameters[] xtis = new Parameters[betai.length];
    	Parameters [][] y = new Parameters[betai.length][iter];    	
    	double [][] acc = new double[betai.length][iter];
    	double [][] loglkd = new double[betai.length][iter];
    	
    
    	Parameters inter;
    	//int lastswap=0;  
    	
    	
    	//boolean stepSelect=true;????????
    	boolean diagnostic=false;
    	
    	boolean[] mixed = new boolean[betai.length];
    	
    	for (int b = 0;b<betai.length;b++){
    		xtis[b] = xin.deepCopy();
    		mixed[b]=false;
    	}


    	boolean n1mult = false;
    	int nswap=0;
    	
    	//tracking of swapping rate per chain
    	double[] ntotswap= new double[betai.length];
    	double[] naccs=new double[betai.length];
    	
    	
      	//logging
    	Writer w1 = log1(betai.length,xtis[0].getParamSize(), dirName);
    	Writer w2 = log2(betai.length,  trialsOptStep, dirName);
    	
    	
    	//Main loop on iteration
    	for (int i = 0;i<iter;i++){	

    		if ((i%100)==0)LOGGER.info("iter "+i);
    		
    		//step tuning on?
    		if (/*stepSelect &&*/ optStep ) {//stepSelect commented out step tuning is stopped when target step is reached ...
    			n1mult = ((i%n1)==0);// done n1 major cycles? if yes, enter step tuning

    			if (n1mult&&(i>0)){
    				
    				double[] rate = new double[betai.length];
    				double[][] oldstep = new double[betai.length][xtis[0].getParamSize()];
    				double[][] step = new double[betai.length][xtis[0].getParamSize()];

    				w1.addVal("iter", i);
    				w1.addVal("swapRate", (nswap*1./n1));
    				
    				//step tuning for all chains
    				for (int b = 0;b<betai.length;b++){   		      				

    					//logger.info("acc.length:"+acc[b].length+" i-n1:"+(i-n1)+" i:"+i);
    					rate[b] = org.apache.commons.math3.stat.StatUtils.mean(acc[b], Math.max(0,(i-n1)), n1);
    					//logger.info("iter "+i+"rate["+b+"]:"+rate[b]);
    					
    		    		w1.addVal(("rate["+b+"]"),rate[b]);
    		    		
    					if (!checkRate(rate[b], targetRate,n1, lambda,  lambdaScale)&& !mixed[b]){
    						oldstep[b] = xtis[b].getSteps();

    						step[b] = stepSelectorX2(diagnostic,betai[b],targetRate,xtis[b], minorCycles, gamma,w2,psim,trialsOptStep, ncyclesOptStep, i, b);//TODO check if we start from xin

    						//System.out.println(FormatUtilities.lineprint("old steps["+b+"]=", oldstep[b], ", "));
    						//System.out.println(FormatUtilities.lineprint("new steps["+b+"]=", step[b], ", "));
    						xtis[b].setSteps(step[b]);	 
    					}  
    					else mixed[b]=true;
    					

    					for (int l=0;l<xtis[0].getParamSize();l++)	w1.addVal(("oldstep["+b+"]["+l+"]"),oldstep[b][l]);
    					for (int l=0;l<xtis[0].getParamSize();l++)	w1.addVal(("step["+b+"]["+l+"]"),step[b][l]);
    				}   



    				w1.vals2File();
    				if (diagnostic) 
    					w1.vals2Term(true);

    				//reset
    				nswap=0;	
    			}  		


    			//stepSelect = false; //stop step 
    		}


    		Object[][] pi = new Object[betai.length][]; 
    		
    		//TODO replace by a call to ChainX ??? with no optimization

    		for (int b = 0;b<betai.length;b++){
    			pi[b] = ptmcmc(xtis[b], betai[b],psim);//
    			xtis[b] = (Parameters)pi[b][0];
    			y[b][i] = xtis[b]; 
    			acc[b][i] = (Double) pi[b][1];
    			loglkd[b][i]=  (Double) pi[b][2];

    		}


    		//swap possible %frequency  when more than 1 chain
    		
    		double u = MCMC.uniformGen(0.,1.)	;
    		if ((u<=swapTargetRate) && (betai.length>2)){
    			//select  betai couple
    			int ibeta = rand.nextInt((betai.length-1));//between 0 and the one before the ultimate index
    			ntotswap[ibeta]++;

    			//accept swap?
    			if (ptAccept(xtis[ibeta], betai[ibeta], xtis[ibeta+1], betai[ibeta+1],psim)){
    				
    				double[] xtisSteps = xtis[ibeta].getSteps();
    				double[] xtis1Steps = xtis[ibeta+1].getSteps();
    				
    				inter = xtis[ibeta].deepCopy();
    				xtis[ibeta] = xtis[ibeta+1].deepCopy();
    				xtis[ibeta+1] = inter.deepCopy();
    				
    				//swap param but keep steps
    				xtis[ibeta].setSteps(xtisSteps);
    				xtis[ibeta+1].setSteps(xtis1Steps);
    				
    				//stepSelect = true;

    				//System.out.print("current rates:");
    				//for (int b=0;b<betai.length;b++)   	System.out.print(", "+StatUtils.mean(acc[b], lastswap, (i-lastswap)));    	
    				//System.out.print(" from "+(i-lastswap)+" samples");	
    				//System.out.println();
    				//logger.info("----------------------------------------swapping @ iter:"+i);

    				//lastswap=i+1;
    				nswap++;
    				naccs[ibeta]++;
    				
    			}
    		}

    	}
    	
    	
    	    	
    	double[] srate=new double[(betai.length-1)];
    	for (int b = 0;b<(betai.length-1);b++) srate[b] = naccs[b]/ntotswap[b];
		String header2 = "swap acceptance rates /chain";
		MCMC.writeAcceptedRates(srate, dirName+"swaprates", header2, "   ");
		
		
    	//return  low temp values ?????
    	this.out = new Object[3*betai.length];
    	for (int b = 0;b<betai.length;b++){
    		this.out[b] = y[b];
    		this.out[(b+betai.length)] = acc[b];
    		this.out[(b+2*betai.length)] = loglkd[b];
    	}

    	w1.flush();
    	w1.close();
    	w2.flush();
    	w2.close();
    	//return out;

   	
    }
    
    /**
	 * MCMC Chain with optional automatic step tuning: can be called alone or from ptChain (from StepSelectorx2)
	 * @param iter number of iteration of the chain
	 * @param xin initial parameter
	 * @param optStep step tuning option
	 * @param betai exponent of likelihood
	 * @param targetRate target acceptance rate
	 * @param psim additional parameters 
	 * @return results
	 * @throws Exception
	 * @throws IOException
	 */
	public static Object[] chainX( Parameters xin, double betai, ParameterSimulation psim, Writer w2, ChainXconfig cfg) throws Exception, IOException{

		
		Parameters xt = xin;
		//Parameters [] y = new Parameters[iter];
		
		double[] steps = xin.getSteps();
		int outiter = cfg.outiter;//number of loops used to tune all the parameters
		int iter = cfg.iter;
		boolean diagnostic = cfg.diagnostic;
	
		//step selection:	
			
		for (int k=0;k<outiter;k++){		
			int iterTest = cfg.iterTest;
			// give a try with smaller chain first to see if step tuning is necessary
			double [] tryacc = (double[]) simpleChainX(iterTest,  xt, betai,psim)[1];
			double rate = StatUtils.mean(tryacc, 0, tryacc.length);
			//logger.info("try rate:"+rate+"");		

			//Set rate for each individual parameter separately 
			int minorCycles = cfg.minorCycles;
			double gamma = cfg.gamma;
			int trialsOptStep = cfg.trialsOptStep;
			int ncyclesOptStep = cfg.ncyclesOptStep;

			//int iterTest = 1;double lambda = 1.;double scale =0.02;
			double lambda = cfg.lambda;
			double scale = cfg.lambdaScale;
			double targetRate = cfg.targetRate;

			if (!checkRate(rate,targetRate,iterTest,lambda,scale))	{
				LOGGER.info("step tuning");
				steps = stepSelectorX2(diagnostic,betai,targetRate,xin,minorCycles, gamma, w2, psim, trialsOptStep, ncyclesOptStep, k, 0);
				LOGGER.info("step tuning done");
			}
			//System.out.println(FormatUtilities.lineprint("new steps", steps, ", "));

			//set new step values
			xt.setSteps(steps);	 
		}			
		
		
		//iterated call to basic mcmc (checked that we are ok here with the copies)
		Object[] res  = simpleChainX(iter,  xt, betai,psim);
		
		
		//return new Object[]{y,acc,loglkd};
		return res;
		
	}
	
	/**
	 * Simple Markov Chain doing an iter number of iterations
	 * @param iter : number of iterations
	 * @param xin
	 * @param betai
	 * @param psim
	 * @return
	 * @throws Exception
	 * @throws IOException
	 */
	public static Object[] simpleChainX(int iter, Parameters xin, double betai, ParameterSimulation psim) throws Exception, IOException{
		Parameters xt = xin;
		Parameters [] y = new Parameters[iter];
		double [] acc = new double[iter];
		double [] loglkd = new double[iter];
		Object[] mcmc;
		
		
		
		//iterated call to basic mcmc (checked that we are ok here with the copies)
		for (int i = 0;i<iter;i++){	 
			mcmc = ptmcmc( xt,betai,psim);
			xt = (Parameters) mcmc[0];
			y[i] = xt; 
			acc[i] = (Double) mcmc[1];
			loglkd[i] = (Double) mcmc[2];
		}
		return new Object[]{y,acc,loglkd};
	}
	
	/**
	     * Core of MCMC. For a parameter, it generates a new pick 
	     * Steps must be set beforehand!!!
	     * @param xt input parameter
	     * @param betai tempering 
	     * @return new parameter with acceptance status (accepted or rejected [old value reused in this case])
	  * @throws IOException 
	  * @throws Exception 
	     */
	    public static Object[] ptmcmc(Parameters xt,double betai, ParameterSimulation psim) throws Exception, IOException{
	    	
	    	
	 	   double acc = 0;
	 	   //proposition for all parameters
	 	   Parameters y;
	 	   
	 	    	   
	 	   Proposal prop = psim.getProposalClass();// new Proposal1beta();
	 	   
	 	   y =  prop.getQProp(xt) ;
	 	   //LOGGER.info("y="+y.getParameter(0));
	 	   /*if (psim==null)  y = lqProp(xt) ;
	 	   else  y = qProp(xt) ;*/
	 		   
	 	   // double a1 = ptPosterior(y, betai, psim)/ptPosterior(xt, betai, psim);
	 	   
	 	   //TEST
	 	   /*for (int p=0;p<xt.nparam;p++){
	 		   if (xt.getParameter(p)<0.) {
	 			   logger.info("erreur signe1");
	 		   }
	 		  if (y.getParameter(p)<0.) {
	 			  logger.info("erreur signe2");
	 		  }

	 	   }*/
	 	   
	 	   
	
	 	   Prior myprior = psim.getPriorClass();
	 	   double priy = myprior.getPrior(y,psim); 
	 	  // double prix = myprior.getPrior(xt,psim);could not be zero or problem with starting values
	 	   
	 	   double r =0;
	 	   Likelihood lkd = psim.getLikelihoodClass();
	 	   double logLkdx = lkd.getLogLikelyhood(xt, psim); 
	 	   double logLkdy =0.;
	 	   //LOGGER.info("logLkdx="+logLkdx);//
	 	   

	 	   //
	 	   if (priy !=0){
	 		   
		 	   logLkdy = lkd.getLogLikelyhood(y, psim);		 	   
	 		   
	 		   double prix = myprior.getPrior(xt,psim);
	 		   double a1 = priy/prix * Math.exp(   betai *  ( logLkdy - logLkdx )    );
	
	 		   // double a2= Qprop(xt)/Qprop(y) ;=1 wwith Gaussian as exp[-(xt-y)^2/sig^2]=exp[-(y-xt)^2/sig^2]
	 		   r = a1;//*a2;
	 		   //LOGGER.info("r="+r);//
	 	   }  
	
	
	
	 	   double u = uniformGen(0.,1.)	;
	
	 	   if ((u<=r) && (r>0)) {//don't return logLkdy if r=0
	 		   //with probab a
	 		   acc=1.;
	 		   return new Object []{y,acc,logLkdy};//add prix x loglkd 
	 	   }
	 	   else    return new Object []{xt.deepCopy(),acc,logLkdx };  //with probab 1-a //add prix x loglkd
	 	  
	 	}



	/** 
	 * check that the targeted rate has been reached
	 * @param rate current rate
	 * @param target rate
	 * @param n1
	 * @param lambda
	 * @param scale
	 * @return
	 */
	public static boolean checkRate(double rate, double target,int n1, double lambda, double scale){
		double epsilon = scale* Math.sqrt(n1*lambda); //at 4% accuracy total
		if (n1*Math.abs(rate-target) <=epsilon ) return true;
		else return false;
		
	}


	/**
	 * tune the step to reach the target acceptance rate
	 * @param diagnostic
	 * @param betai
	 * @param targetRate
	 * @param xin
	 * @param minorCycles
	 * @param gamma
	 * @param w2
	 * @param psim
	 * @param trials
	 * @param ncycles
	 * @param iterNum
	 * @param chainNum
	 * @return
	 * @throws Exception
	 * @throws IOException
	 */
	public static double[] stepSelectorX2(boolean diagnostic, double betai, double targetRate, Parameters xin, int minorCycles, double gamma, Writer w2,ParameterSimulation psim, int trials, int ncycles, int iterNum, int chainNum) throws Exception, IOException  {
		int nparam=xin.getParamSize();
		
		//int trials = 11;//11 number of chains to run : 1 per trial step
		double [] optistep = new double[nparam] ;//array of optimal steps: 1 per parameter
		Parameters xt;
		//int ncycles =2;//1;//20 //could use ncycles to get a sample of steps for each parameter but too much resources needed...
	 	
		//step selection for each varying parameter independantly 
		for (int p =0;p<nparam;p++){
			if (xin.getVarying(p)){

				xt = xin.deepCopy();   

				double optis2=xt.getStep(p);  

				for (int q =0; q<ncycles; q++){//only 1 q    	   		

					double[] acc = null;    		

					double[] accr =  new double[trials ];
					double[] sigr =  new double[trials ];

					List<Double> mr= new ArrayList<Double>() ;
					List<Double> ms= new ArrayList<Double>() ;
					//Parameters[] yl;
					

					int stp2small = 0;
					int stp2big =0;
					double slope = -100.;
					double offset = 0.;


					//loop where trying all the steps for 1 param
					// runs a number of "trials" of chains 
					for (int i=0;i<trials;i++) {

						//running 1 chain to calculate rate ...
						xt = xin.deepCopy();        					
						sigr[i] = optis2*10./Math.pow(2., i);   //10     			

						xt.setSteps(p, sigr[i]);//very important here: for each trial, a new  step
						
                        //perform an minorCycles number of iterations
						Object[] ch =  simpleChainX(minorCycles, xt, betai, psim);

						//yl = (Parameters []) ch[0];
						//yp= Parameters.extractor(yl,p);//?
						acc = (double[]) ch[1];        			

						//calculating rate
						accr[i]= StatUtils.mean(acc.clone(), 0, acc.length);

						//need to filter out 0 and 1 for logit function	
						if (accr[i]>0.){
							if (accr[i]<1.) {			
								//store step and rate
								mr.add(accr[i]);
								ms.add(sigr[i]);
							}else stp2small++;//if step =100%    , step too small   

						}
						else stp2big++;//if rate =0 ... step is too big	

					}

					//	System.out.println(FormatUtilities.lineprint("sigma=", sigr, ", "));        		
					//	System.out.println(FormatUtilities.lineprint(" rate=", accr, ", "));   


					w2.addVal("iter", iterNum);
					w2.addVal("chain", chainNum);
					w2.addVal("param",p);
					w2.addVal("cycle",q);

					for (int l=0;l<trials;l++) w2.addVal(("step["+l+"]"),sigr[l]);
					for (int l=0;l<trials;l++) w2.addVal(("rate["+l+"]"),accr[l]);



					String mode="o";
					int stepLim = 5; //limit on the number of too big or too small steps
					if ((stp2big<stepLim)&&(stp2small<stepLim)) {//if not too many small steps, can calculate a optimal step value from data points
						double[] sol =   stepCalculator2(targetRate, mr.toArray(), ms.toArray(),diagnostic,p);
						optis2 = sol[0];
						slope = sol[1];
						offset = sol[2];

						/*
	       		 double[] resint = new double[3];
	       		 resint = stepCalculator3(targetRate, mr.toArray(), ms.toArray(),diagnostic,offset, slope, trials);
	       		double optis3 = resint[0];
	       		offset = resint[1];
	       		slope = resint[2];*/
					}
					if (stp2big>=stepLim) {//if too many big steps, must divide step by 10
						optis2 =optis2/10.;//logger.info("optimal step divided");
						mode ="d";//for info
					}
					if ((stp2small>=stepLim)||(slope> -0.5  )) {// too many small step or slope >-0.5
						//when more than 1 param, rate can saturate under  100% even for tiny step
						//in this case, slope coefficient get closer to 0 from negative side
						optis2 =optis2*10.;//logger.info("optimal step multiplied");
						mode="m";//for info
					}

					//logger.info("optimal step for param:"+p+"="+optis2+" @cyle:"+q);


					w2.addVal("opti", optis2);
					w2.addVal("mode", mode); 	
					w2.addVal("slope", slope);
					w2.addVal("offset", offset); 	

					w2.addVal("type", "step");            	
					if (diagnostic) w2.vals2Term(true, "rate");
					w2.vals2File("rate");

					w2.addVal("type", "rate");    
					if (diagnostic) w2.vals2Term(true, "step");
					w2.vals2File("step");
				}
				//maximum damping to avoid to big changes
				if (gamma!=0.) {

					if (xin.getStep(p)<optis2) optistep[p] = Math.min(gamma* xin.getStep(p),optis2 );
					else optistep[p] = Math.max(1./gamma* xin.getStep(p),optis2 );

				}
				//if want to keep last param value
				//optistep[p] = optis2;


			}
		}
		return optistep;	
	
	}


	/**
	 * calculate the optimal step to reach the target acceptance rate from interpolation  (tested)
	 * @param targetRate target acceptance rate
	 * @param rate array of rate obtained from chains with different step
	 * @param step array of proposed step 
	 * @param diagnostic
	 * @param p
	 * @return
	 * @throws Exception
	 * @throws IOException
	 */
	public static double[] stepCalculator2(double targetRate, Object[] rate, Object[] step, boolean diagnostic, int p ) throws Exception, IOException{
		
		double[] logitrate = new double[rate.length];
		double[] lstep = new double[rate.length];
		double[] mrate = new double[rate.length];
		double[] mstep = new double[rate.length];
		
		
		//regression: logit(rate) = slope * log(rate) + offset
		for (int i=0;i<rate.length;i++ ){
			
			mrate[i]= (Double) rate[i];
			mstep[i]=(Double)step[i];
			
			logitrate[i] = logit(mrate[i]);
			if (mstep[i]>0.) lstep[i] = Math.log(mstep[i]);
			else throw new IllegalArgumentException("A step value must be strictly larger than 0");//assume step != 0 which should be the case
		}
		
		
		//linear interpolation: the choice of the (x,y) vs (y,x) is important
	
		/***old method:
		final LinearRegression.Result result = LinearRegression.compute(lstep,logitrate);
		double offset = result.getOrigin();
		double slope = result.getSlope();***/		
		
		SimpleRegression simpleRegression = new SimpleRegression(true);

		for (int i=0;i<rate.length;i++ ){
			simpleRegression.addData(lstep[i],logitrate[i]);
		}
		// querying for model parameters
		double offset =  simpleRegression.getIntercept();
		double slope = simpleRegression.getSlope();

		
		//optimal step from regression
		double optiStep = Math.exp((logit(targetRate)-offset)/slope  ); 
		
		//Possible improvement:  check slope error or check if target rate in regression data range.
		
		//check what is done
		if (diagnostic) {
			LOGGER.info("offset2:"+offset+" slope2:"+slope);
			LOGGER.info("logit(targetRate):"+logit(targetRate));
	
	
			//control: calculate rates from regression law
			double[] lcalcRate = new double [rate.length];
			double[] calcRate = new double [rate.length];
			for (int i=0;i<rate.length;i++ ){
				lcalcRate[i]= ((slope*lstep[i]+offset) );
				calcRate[i]= Math.exp(lcalcRate[i]) / ( Math.exp(lcalcRate[i])+1.);
	
			}
		}
		
		return new double[]{optiStep,slope,offset};
	}


	/**
	 * Test if condition for parameter chain swapping is ok
	 * @param xti
	 * @param betai
	 * @param xti1
	 * @param betai1
	 * @param psim
	 * @return
	 * @throws Exception
	 * @throws IOException
	 */
	public static boolean ptAccept(Parameters xti,double betai,Parameters xti1,double betai1 , ParameterSimulation psim) throws Exception, IOException{
		
		//prior compensates each other in ratio but need to calculate them for the case equal to 0
		Prior myprior = psim.getPriorClass();
		double pri3 = myprior.getPrior(xti,psim); 
		double pri4 = myprior.getPrior(xti1,psim); 
	
		double r =0;
		if ((pri3*pri4) !=0){//division by zero but null would mean posterior is 0 (or r=0)
	
			//NOTE pri1 = prior(xti1,psim)=pri4 and  pri2 = prior(xti,psim)=pri3
	
			Likelihood lkd =  psim.getLikelihoodClass();
			double a1 = Math.exp(   betai *  (lkd.getLogLikelyhood(xti1, psim) - lkd.getLogLikelyhood(xti, psim) )    );//times ratio pri1/pri3 
			double a2 = Math.exp(   betai1 * (lkd.getLogLikelyhood(xti, psim) - lkd.getLogLikelyhood(xti1, psim) )    );//*pri2/pri4 
			r = a1*a2;
	
		}
	
		double u = uniformGen(0.,1.)	;
	
		if (u<=r)  {
			//swap
			return true;
		}
		return false;
	 	
	 }


	/**
	 * logit function (tested)
	 * @param x
	 * @return
	 */
	public static double logit(double x){
		
		return Math.log(x/(1.-x));
	}


	/**
	 * range limited & shifted uniform random generator
	 * @param mu offset
	 * @param sigma range
	 * @return
	 */
	public static double uniformGen(double mu,double sigma){
		return mu + sigma* rand.nextDouble();
	}


	/**
	 * Gaussian distributed sample generation
	 * @param mu
	 * @param sigma
	 * @return
	 */
	public static double gaussGen(double mu,double sigma){
		return mu + sigma* rand.nextGaussian();
	}


	/**
     * Calculate each chain averaged acceptance over the number of iterations
     * @return average rate[chain]
     */
    public double[] getChainAcceptanceRates(){
    	double[] acc;
    	double[] rate = new double[lbetai.length];
    	
		for (int nch=0;nch<lbetai.length;nch++){
			
			acc = (double[]) this.out[(nch+lbetai.length)];				

			rate[nch] = StatUtils.mean(acc, 0, acc.length);
		}
    	return rate;
    }
    
    /**
     * Getter for the log likelihood values for each chain 
     * @return log likelihood[chain][iteration]
     */
    public double[][] getChainLogLikelihoods(){
    	double[][] loglkd =  new double[lbetai.length][];
		for (int nch=0;nch<lbetai.length;nch++)		loglkd[nch]  = (double[]) this.out[(nch+2*lbetai.length)];			
    	return 	loglkd;
    }
    
    
    /**
     * Return array of samples of parameters for each iteration and each chain.
     * The list of parameter value is given for each chain
     * @return sample value [iteration][parameter order]
     */
    /* not clear if in use
    public double[][] getAllChainSamples(){
    	double[][] res =null;
    	double[][] finres =null;
    	Parameters[] yp;
    	
		for (int nch=0;nch<lbetai.length;nch++){
			yp = (Parameters[]) this.out[nch];
			
			
			res = MCMC.extractor(yp);
			if (nch==0) finres= Array2DUtilities.arrayClone2D(res);			//???
			else 				finres = Array2DUtilities.arrayMergeCol(finres,res);
			
		}
    	return finres;
    }   
    */
    /**
     * Return array of samples of parameters for each iteration for 
     * a specific chain.
     *
     * @param chain
     * @return array of parameters
     */
    public Parameters[] getChainParameters(int chain){

    	Parameters[] yp;    	

    	yp = (Parameters[]) this.out[chain];			

    	return yp;
    }   
   
    /**
     * convert a list of parameters into a 2D array
     * one lineper iteration* nparameters column for samples + nparameters columns for steps
     * @param yl
     * @return 2D array of parameters
     */
    public static double[][] extractor(Parameters[] yl){
    	int iter  = yl.length;
    	int nparam = yl[0].getParamSize();
    	double[][] param = new double[iter][2*nparam];
    	
    	for (int i=0;i<iter;i++)  {

    		for (int j=0;j<nparam;j++){
    			param[i][j] = yl[i].getParameter(j); 
    			param[i][j+nparam] = yl[i].getStep(j);
    		}
    	}
    	return param;
    } 
    
 
	public static Parameters[] getChainXParameters(Object[] data){
    	Parameters[] p = (Parameters[]) data[0];
    	return p;
    }
    public static double[] getChainXAcceptance(Object[] data){
    	double[] a = (double[])data[1];
    	return a;
    }
    public static double[] getChainXlogLkd(Object[] data){
    	double[] l = (double[])data[2];
    	return l;
    }
    
    /**
     * sample median for all parameters
     * @param endBurn
     * @param chain
     * @return array of median values
     */
	public double[] getSamplesMedian(int endBurn, int chain) {
		Parameters[] yp;		
		
		yp = getChainParameters(chain);
    	int psize =yp[0].getParamSize();
    	
		// first flavour of results
    	double[] median = new double[psize];
    	
    	LOGGER.info("for median, end of burning period set to ="+endBurn+" iterations");

		for(int n =0;n<psize;n++){
			double[] y =new double[yp.length] ;
			y = Parameters.extractor(yp,n);						
					
			DescriptiveStatistics stats = new DescriptiveStatistics();
			
			for (int i = endBurn; i < y.length ;i++){
				stats.addValue(y[i]);
			}

			//double mean = stats.getMean();		
			
			LOGGER.info("median["+n+"]:"+median[n]);		

		}
		
		
		return median;
	}
	   public static double getSamplesMap(int chain, double xmin, double xmax ,double[] samples, double [] loglkd){
		double lkdmax = -1.*Double.MAX_VALUE ;
		int index = 0;
		for(int n =0;n<loglkd.length;n++){
			if (   (loglkd[n]>lkdmax) && ( (samples[n]>xmin) && (samples[n]<xmax)   )  ) {
				lkdmax = loglkd[n];
				index = n;
			}
			
		}
		return samples[index];
	}


	/**
     * Getter for samples  for all parameters for a chain
     * @param endBurn
     * @param chain
     * @return array of median values
     */
	public double[] getSamples( int chain, int param) {
		Parameters[] yp;		
		
		yp = getChainParameters(chain);

		double[] y = Parameters.extractor(yp,param);						
		return y;
	}
    /*NOT READY
	public void samplePlot(int chain, int n, String filePathName, String fileName) throws Exception {
		
		
		double[] y ;
		y = getSamples(chain, n);				

		
		double[] x = Idl.dindgen(y.length);

		PlotXYData plot = new PlotXYData(("mcmc:"+n),"iteration" , "Log",
				AxisType.LINEAR, AxisType.LOGARITHMIC);
		plot.addData(x, y);

		if (filePathName == null ) plot.display((short) 700);
		else plot.save(filePathName, fileName, RenderTarget.PNG, (short)700);	
		
	}
	
	public void sampleHistPlot(int chain, int n, String filePathName, String fileName) throws Exception {
		double[] y ;
		y = getSamples(chain, n);				
		
		
		
		int nrBins = 50;
		PlotHist ploth = new PlotHist( nrBins, ("mcmc:"+n) );
		ploth.addData(y, false, ("mcmc:"+n), Color.BLUE);		
		if (filePathName == null ) ploth.display((short) 700);
		else ploth.save(filePathName, fileName, RenderTarget.PNG, (short)700);	
	}

    */
	/**
	 * write MCMC proposals to a file, one column for each parameter
	 * @param data results
	 * @param ASCIIfileNameExtensless
	 * @param header
	 * @param delimiter
	 */
	public static void writeMcmcData(double[][] data,String ASCIIfileNameExtensless, String header, String delimiter){
		String s = ASCIIfileNameExtensless
		+ ( new SimpleDateFormat ("yyMMdd-HHmmss") ).format(Calendar.getInstance().getTime())+".dat";
		
		
			File targetFile = new File(s);
	
			try {
				PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(targetFile, true)));
				pw.printf("%s \n",header);
				
				for (int i = 0; i < data.length; i++) {
					for (int j = 0; j < data[0].length; j++) {
	
						pw.printf("%14.6e %s", data[i][j], delimiter);
					}
					pw.println();
	
				}
				
				pw.close();
	
			} catch (FileNotFoundException fe) {
				LOGGER.severe((" error in writeMcmcData:"+ fe.getCause()));
				fe.printStackTrace();
				System.exit(0);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		
			
	}

	public static void writeTestLikeData(double[][] data,String ASCIIfileNameExtensless,  String delimiter){
		String s = ASCIIfileNameExtensless
		+ ( new SimpleDateFormat ("yyMMdd-HHmmss") ).format(Calendar.getInstance().getTime())+".dat";
		//double[][] datatest = Array2DUtilities.mirrorAC(Array2DUtilities.mirrorAL(data)); WHY?
		
			File targetFile = new File(s);
	
			try {
				PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(targetFile, true)));
				
				
				for (int i = 0; i < data.length; i++) {//was datatest
					for (int j = 0; j < data[0].length; j++) {//was datatest
	
						pw.printf("%14.6e %s", data[i][j], delimiter);
					}
					pw.println();
	
				}
				
				pw.close();
	
			} catch (FileNotFoundException fe) {
				LOGGER.severe((" error in writeMcmcData:"+ fe.getCause()));
				fe.printStackTrace();
				System.exit(0);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		
			
	}

	public static void writeAcceptedRates	(double[] data,String ASCIIfileNameExtensless, String header, String delimiter){
		String s = ASCIIfileNameExtensless
		+ ( new SimpleDateFormat ("yyMMdd-HHmmss") ).format(Calendar.getInstance().getTime())+".dat";
		
		
			File targetFile = new File(s);
	
			try {
				PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(targetFile, true)));
				pw.printf("%s \n",header);
				
				for (int i = 0; i < data.length; i++) pw.printf("%e %s", data[i], delimiter);					
			
				
				pw.close();
	
			} catch (FileNotFoundException fe) {
				LOGGER.severe((" error in writeAcceptedRates:"+ fe.getCause()));
				fe.printStackTrace();
				System.exit(0);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		
	
	}
}
