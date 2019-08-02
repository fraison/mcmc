package algo.mcmc;



import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Locale;
import java.util.ResourceBundle;


/**
 * This class provides parameter for ChainX especially for call to the step selector (stepSelectorX2)
 * @author fraison
 *
 */

public class ChainXconfig {
	
	//specific to ptChainX
	 double swapTargetRate;
	 public boolean optStep;
	 public int majorCycles;
	 public String dirName;
	 public double[] betai; 
	
	//Shared ptChainX & ChainX -> READ from properties in constructor
 	 public int iter;
	 public int minorCycles;//= 100;//100 50	for StepSelectorX2
	 public double gamma ;//= 30; for StepSelectorX2
	 public int trialsOptStep ;//= 11;//11 number of chains to run : 1 per trial step;for StepSelectorX2
	 public int ncyclesOptStep ;// = 2;//1;//20 //could use ncycles to get a sample of steps for each parameter but too much resources needed...;for StepSelectorX2
	 public double lambda ;//= 0.25;for checkRate
	 public double lambdaScale ;//= 1.5;for checkRate
	 public double targetRate;//for checkRate
	 public boolean diagnostic;//for plots and xterm logging of step tuning
	 
	 //specific to ChainX -> Not READ from properties
	 public int iterTest = 50 ;////default value : number of iteration in assessment of acceptance rate prior to tuning NOT READ from properties
	 public int outiter = 1 ;////default value : number of loops used to tune all the parameters

	



	 /**
	  * constructor initiating values from property file
	  * @param props
	  * @throws GaiaConfigurationException
	  */
	 public ChainXconfig(String props) throws Exception{
		 //PropertyLoader.load(props);
		 Locale locale = new Locale("en", "US");

		 ResourceBundle resource = ResourceBundle.getBundle(props, locale);
		 
	     //specific to ptChainX
		 this.swapTargetRate = Double.parseDouble(resource.getString("mcmc.algo.swapTargetRate"));
		 this.optStep = Boolean.parseBoolean(resource.getString("mcmc.algo.optimizeStep"));
		 this.majorCycles = Integer.parseInt(resource.getString("mcmc.algo.majorCycles"));
		 this.dirName = resource.getString("mcmc.algo.dirName");
		 String[] betaString = resource.getString("mcmc.algo.betai").split(",");
		 //this.betai = Arrays.stream(betaString).mapToDouble(Double::parseDouble).toArray(); java 8 
		 this.betai = new double[betaString.length];
		 for (int i=0;i<betaString.length;i++){ this.betai[i]= Double.parseDouble(betaString[i]);}
				 
		 
		 
		 
		//Shared with ptChainX & ChainX
		 this.iter = Integer.parseInt(resource.getString("mcmc.algo.iter"));    	 
		 this.minorCycles = Integer.parseInt(resource.getString("mcmc.algo.minorCycles"));
		 this.gamma = Double.parseDouble(resource.getString("mcmc.algo.gamma"));	
		 this.trialsOptStep = Integer.parseInt(resource.getString("mcmc.algo.ntrialsOptStep"));
		 this.ncyclesOptStep = Integer.parseInt(resource.getString("mcmc.algo.ncyclesOptStep"));
		 this.lambda = Double.parseDouble(resource.getString("mcmc.algo.lambda"));
		 this.lambdaScale = Double.parseDouble(resource.getString("mcmc.algo.lambdaScale"));  
		 this.targetRate = Double.parseDouble(resource.getString("mcmc.algo.targetRate"));
		 this.diagnostic = Boolean.parseBoolean(resource.getString("mcmc.algo.diagnostic"));
		 
    	
	 }
	 
	 /**
	  * save run configuration
	  * @throws IOException
	  */
	 public void saveXconfig(String path) throws IOException{
		 String s;
		 if (path==null) s = this.dirName+"config_" + ( new SimpleDateFormat ("yyMMdd-HHmmss") ).format(Calendar.getInstance().getTime())+".dat";
		 else s = path+"config_" + ( new SimpleDateFormat ("yyMMdd-HHmmss") ).format(Calendar.getInstance().getTime())+".dat";

		 File targetFile = new File(s);


		 PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(targetFile, true)));
		 //pw.printf("%s \n",header);

		 pw.printf("%s %f\n", "swapTargetRate = ",this.swapTargetRate);
		 pw.printf("%s %s\n","optStep = ",this.optStep);
		 pw.printf("%s %d\n","majorCycles = ",this.majorCycles);
		 pw.printf("%s %s\n","dirName = ",this.dirName);
		 pw.printf("%s ","betai ="); for (int i = 0; i <betai.length;i++) pw.printf("%e %s",this.betai[i]," ,");pw.println();


		 pw.printf("%s %d\n","iter = ",this.iter);
		 pw.printf("%s %d\n","minorCycles = ",this.minorCycles);
		 pw.printf("%s %f\n", "gamma = ",this.gamma);
		 pw.printf("%s %d\n","trialsOptStep = ",this.trialsOptStep) ;
		 pw.printf("%s %d\n","ncyclesOptStep = ",this.ncyclesOptStep) ;
		 pw.printf("%s %f\n", "lambda = ",this.lambda) ;
		 pw.printf("%s %f\n", "lambdaScale = ",this.lambdaScale) ;
		 pw.printf("%s %f\n", "targetRate = ",this.targetRate) ;
		 pw.printf("%s %s\n","diagnostic = ",this.diagnostic);


		 pw.printf("%s %d\n","iterTest = ",this.iterTest);
		 pw.printf("%s %d\n","outiter = ",this.outiter);			


		 pw.close();


		
	 }
}
