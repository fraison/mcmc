# MCMC with parallel tempering and automatic step tuning in Java

This code is corresponding to the descritption on my [blog](https://multipass193830787.wordpress.com/mcmc-with-parallel-tempering-and-automatic-step-tuning-in-java/).

## Implementation

The MCMC class implements the algorithms to perform parameter inference with Markov Chain Monte Carlo with parallel tempering algorithm.

The class *ChainXconfig* is a convenience class used to collect and store the MCMC control parameters specified in "mcmc.properties".

The abstract class *Parameters* has 3 members that need to be set in the constructor:

- the initial values of the model parameters to be inferred stored in a one dimension array of "double".
- the initial step values of the model parameters to be inferred stored in a one dimension array of "double". This corresponds to the standard deviations of the proposal distribution and the order of values should correspond to the array of initial values above.
- the parameters to be updated through an array of "boolean". If a parameter should be updated, then its value should be "true". Again, the values in this array should be in the same order of the two previous arrays.

The interface *Proposal* corresponds to the proposal distribution q(Y|X) in equation 2 in my blog. It declares a method *getQProp(Parameters c)* which is supposed to return the new proposed Parameters from the previous Parameters object. The method gaussGen(double mu, double sigma) can be used to implement a multivariate normal Gaussian. See for example the class *ProposalAnybeta* in the package *mcmc.progs*.

The interface *likelihood* corresponds to the likelihood distribution in equation 2 in my blog. It declares a method *getLogLikelyhood(Parameters x, ParameterSimulation psim)* that returns a "double".

The interface *Prior* corresponds the prior distribution p(X|I) in equation 1. It declares a method *getPrior(Parameters x, ParameterSimulation psim)* that returns a "double".

The abstract class *ParametersSimulation* is a convenience class that carries information necessary for the likelihood calculation: data set to which simulated data from the model are compared, additional parameters to make simulated data from the model, etc.

The class MCMC contains the principal methods to manage the chains like "ptChainX" and the algorithms like the automatic step tuning ("stepSelectorX2").


## Properties for MCMC control parameters

- mcmc.algo.swapTargetRate is the target value for chain swapping acceptance rate. Default value : 0.0333; this means once every 30 iterations, a random value between [0,1] is generated and if its value is less then "swapTargetRate", we enter a loop to start a potential parameter swapping between 2 chains. 
- mcmc.algo.optimizeStep if "true", the chain will perform the automatic step tuning. Otherwise, step will remain equal to their initial value. Default: true.
- mcmc.algo.majorCycles number of iterations after which a step tuning is started. Default value: 1500.
- mcmc.algo.dirName is the full directory name for output in the workspace under the project root. Default: "data/output/".
- mcmc.algo.betai is the 1D array of exponents *beta_i* of Likelihood in parallel tempering (the number of exponents defines the number of tempering chains). The array of "double" must be in growing order and the last one should be 1 (the lowest temperature in the parallel tempering scheme). Default values: 1. This means a single chain without tempering (standard MCMC). 
Example: 1e-3, 2.5e-3, 4.7e-3, 7.7e-3, 1.5e-2, 2.575e-2, 5.05e-2, 7.525e-2, 1.3e-1, 2.5e-1, 3.5e-1, 5.e-1, 6.5e-1, 1..
- mcmc.algo.iter is the number of iterations. Default value: 100000. Corresponds to the final number of available samples per parameter in the posterior distribution.
- mcmc.algo.minorCycles is the number of minor cycles e.g. the number of iterations of each chain in the step tuning module. Default value: 50.
- mcmc.algo.gamma is the maximum value of damping factor on step selection. Default value: 1.6. 
- mcmc.algo.ntrialsOptStep is the number of step values to calculate the optimal step (1 chain per step).Default value: 11.
- mcmc.algo.ncyclesOptStep is the number of cycles of ntrials to use per parameter to find the optimal step; This is loop doesn't include the damping. Default value: 1.
- mcmc.algo.lambda is a damping factor to calculate the tolerance on the target acceptance rate. Default value: 0.25. 
- mcmc.algo.lambdaScale is a scaling factor to calculate the tolerance on the target acceptance rate. Default value: 1.5.
- mcmc.algo.targetRate is target value for acceptance rate. Default value: 0.20 (warning: NOT a percentage).
- mcmc.algo.diagnostic trigger plotting and xterm logging of step tuning phase. Default: false. Should be used only for debugging as it generates a substantial amount of data.

## Operation

Here is the list of operations to perform to run the algorithm. The list has an order.

1. Set the MCMC control parameters in file "mcmc.properties".
2. Choose the initial parameters (1D array of "double").
3. Choose the initial values of the step for each parameter (1 dimension array of "double").
4. Choose the updated parameters (1D array of "boolean").
5. Instantiate an object "prop" from the implementation of the class *Proposal*.
6. Instantiate an object "lkd" from the implementation of the class *likelihood*.
7. Instantiate an object "prior" from the implementation of the class *Prior*.
8. Instantiate an object "xin" from a child class of the *Parameter* class.
9. Instantiate an object "psim" from a child class of the *ParameterSimulation* class through the constructor *ParameterSimulation(prop, prior, lkd)*.
10. Instantiate an object "cfg" from the *ChainXconfig* class.
11. Instantiate an object from the *MCMC* class.
12. Call the method "ptChainX(xin, psim, cfg)" from previous object.
13. Obtain output from the available getters of the *MCMC* class.

## Output
The samples are the final product of the software. However, a detailed logging has been organized as the tuning of parameters can be complex.

### Data

The output of the "ptChainX" is an array of objects composed of:

- the arrays, for all chains, of arrays of objects *Parameters[iteration]* for all iterations,
- the arrays, for all chains, of arrays giving the status of proposed *Parameters[iteration]* for all iterations,
- the arrays, for all chains, of arrays giving the logarithm of the likelihood of the model for *Parameters[iteration]* for all iterations.

To facilitate the extraction of the data, some getters have been written:

- "getAllChainSamples()": extracts samples and step values in an array of "double". The array stores for each iteration, the full list of parameters and then the full list of steps for each chain. The table can be easily saved into a text file with the method "writeMcmcData" from MCMC.
- "getChainAcceptanceRates()": calculates the average acceptance rate for each chain. This can be saved into a text file with method "writeAcceptedRates"from MCMC.
- "getChainLogLikelihoods()": for each iteration, extracts the logarithm of the likelihood for each chain. This can be saved into a text file with method "reverse" from Array2DUtilities followed by method "writeMcmcData" from MCMC.

Marginal distributions can be visualized through the histogram of each parameter samples and parameters can be simply estimated from the data. For example, the median of marginal distribution can be calculated from the median of the samples for each parameter. Also, the maximum a posteriori (MAP) can be calculated taking the samples values corresponding to the iteration providing the largest value of the log Likelihood. Both methods give good results provided that the samples from the burn-in period are discarded. However these methods can be put into questions when the distribution is not uni modal or when the chains didn't reach an equilibrium. 

### Loggers

The logger are an help to understand the impact of the control parameters on the chain. The file "log1_{date}.dat" gives information related to the global control of the tempering. The file "log2_{date}".dat" gives information related to the control of the step tuning.

## Example

The method "testptChainX" in class *MCMCTest* from the package *mcmc.algo.mcmc.test* can give a good start on how to use the parallel tempering in a simple setting (1 parameter) but in a complex case (multi peaked Gaussian).

