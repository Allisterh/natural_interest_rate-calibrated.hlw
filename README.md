# calibrated.hlw

This code is used to estimate the natural rate of interest (r-star) using a calibrated Holston, Laubach and Williams (2017) multi-stage maximum likelihood estimation procedure, where the likelihood function is computed by the Kalman filter. For more details on the state-space model, refer to their paper.

The function run.hlw.estimation.R is called by run.hlw.R once for each economy. It takes as inputs the key variables for the given economy: log output, inflation, and the real and nominal short-term interest rates, as well as the specified constraints on ar and by. It calls the programs rstar.stageX.R to run the three stages of the HLW estimation. Additionally, it calls the programs median.unbiased.estimator.stageX.R to obtain the signal-to-noise ratios λg and λz.

The programs unpack.parameters.stageX.R set up coefficient matrices for the corresponding state-space models for the given parameter vectors. In all stages, the constraint on b ≥ 0.025 is imposed. In stages 2 and 3, a ≤ −0.0025. These constraints are labeled as a3.constraint and b2.constraint, respectively, in the code.

For more detail on each stage of the estimation procedure, consult HLW (2017). The programs rstar.stageX.R run the model. The function median.unbiased.estimator.stage1.R computes the exponential Wald statistic of Andrews and Ploberger (1994) for a structural break with unknown break date from the first difference of the preliminary estimate of the natural rate of output from the stage 1 model to obtain the median unbiased estimate of λg.

The function median.unbiased.estimator.stage2.R applies the exponential Wald test for an intercept shift in the IS equation at an unknown date to obtain the median unbiased estimate of λz, taking as input estimates from the stage 2 model.

Within the program kalman.states.R, the function kalman.states() calls kalman.states.filtered() and kalman.states.smoothed() to apply the Kalman filter and smoother. It takes as input the coefficient matrices for the given state-space model as well as the conditional expectation and covariance matrix of the initial state, xi.tm1tm1 and P.tm1tm1, respectively. kalman.states.wrapper.R is a wrapper function for kalman.states.R that specifies inputs based on the estimation stage.

The function kalman.log.likelihood.R takes as input the coefficient matrices of the given state-space model and the conditional expectation and covariance matrix of the initial state and returns the log likelihood value and a vector with the log likelihood at each time t. log.likelihood.wrapper.R is a wrapper function for kalman.log.likelihood.R that specifies inputs based on the estimation stage.

The function kalman.standard.errors.R computes confidence intervals and corresponding standard errors for the estimates of the states using Hamilton’s (1986) Monte Carlo procedure that accounts for both filter and parameter uncertainty. The function calculate.covariance.R calculates the covariance matrix of the initial state from the gradients of the likelihood function. The function format.output.R generates a dataframe to be written to a CSV containing one-sided estimates, parameter values, standard errors, and other statistics of interest.

For each economy, the final section of run.hlw.R reads in prepared data, runs the HLW estimation by calling run.hlw.estimation.R, and saves one- sided estimates and a spreadsheet of output.
