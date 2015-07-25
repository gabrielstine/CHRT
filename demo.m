%
% Examples of using Choice & Reaction Time (CHRT) Matlab toolbox to fit and 
% generate 1D Choice & Recation Time data.
%
% By Shadlen Lab, Columbia University in the city of New York.
%
%
% In this Choice & Reaction Time toolbox, 1D fitting data is sorted by trial 
% and listed by row. For each trial, it is arranged as
%
%   [signed coherence, choice, reation time]. 
%
% A sample data has been provided for fitting. First load and check the sample 
% data. Place a break point at line blow and step through this demo.

load data.mat
whos
data

% This sample data is a 3-column matrix. The first column is 'signed
% coherence' in the range of [-0.512, 0.512], of which negative value means
% the dots moves toward left and positive value means the dots moves toward
% right. The second column is 'choice', of which '1' means choosing right
% while '0' means left. The third column is 'reaction time' in unit of
% second.

%% Logit Regression Fit.
% For the provided sample data, the unique signed coherence and the 
% corresponding probability of making rightward choice satisfy the Logit
% Regression profile,
%
%   log(pr/(1 - pr)) = Xb = b1 + b2*x1 = b2 * (b1/b2 + x1)
%   x1 = unique signed coherence,
%   pr = probability of making rightward choice.
%
% To apply Logit Regression Fit, simply call LOGITFIT function with 
% sample data as input. This function will calculate the probability of
% making rightward choice and call GLMFIT function to perform the actual fit.

[b,dev,stats,prData] = logitFit(data)

% LOGITFIT function has 4 output variables, of which 
%   b = a vector of coefficient estimates, 
%   dev = the deviance of the fit,
%   stats = the structure containing fitting fields,
%   prData = the calculated proportion rightward choice data.
%
% Of those 4 output variables, b, dev & stats are the same as the output
% variables of GLMFIT function. In the above Logit Regression profile, b1/b2 
% is called 'logit derived coherence bias' and can be used in the following 
% computations. 

logitDerivedCohBias = b(1) / b(2)

% prData can be used in plot. To visualize Logit Regression fit curve and 
% the proportion rightward choice data, call the Logit Regression Fit plot
% function by supplying the proportion rightward choice data and the
% coefficient estimate vector. In the plot, the blue dots stand for the
% acutal data supplied and the cyan color curve for the fitting data. The 
% plot function also returns the figure handle which can be used to modify
% the figure. User should feel free to generate figures with their own
% preferences.

fh1 = logitFitPlot(prData, b)

%% Weibull fit.
% The unsigned coherence and the corresponding probability of making correct 
% choice also satisfy Weibull profile,
%
%   pc = 1 - 0.5 * exp(-(x/alpha).^beta)   
%
% alpha & beta are the threshold and slope parameters.
% To apply Weibull fit, call WEIBULLFIT function by supplying 'data' as
% the input. This WEIBULLFIT function returns 5 outputs which are
%   alpha = threshold parameter,
%   beta = slope parameter,
%   llik = log likelihood of obtaining data given the fit,
%   abse = 2-vector with standard errors for alpha and beta,
%   pcData = proportional correct choice data of which the 3 columns are
%       unsigned coherence, probability of correct choice and number of
%       total trials.

[alpha,beta,llik,abse,pcData] = weibullFit(data)

% The Weibull fit data can be visualized by calling the designed plot
% function WEIBULLFITPLOT with the threshold parameter 'alpha', the slope
% parameter 'beta' and the condensed proportional correct choice data
% 'pcData'. In the plot, the blue dots stand for the acutal data supplied 
% and the cyan color curve for the fitting data generated from Weibull
% profile. This plot function also returns a function handle for plot
% editing. 

fh2 = weibullFitPlot(alpha,beta,pcData)

%% Flat Bound fit.
% As described in book "Bayesian Brain: Probabilistic Approaches to Neural 
% Coding" edited by Kenji Doya, Shin Ishii, Alexandre Pouget and Rajesh PN.,
% Rao, Chapter 10, "The speed and accuracy of a simple perceptual decision:
% A mathematical primer", the probability of choosing the
% positive/rightward direction and the decision time can be expressed as a 
% function of the motion strength (coherence),
%
%   p = 1 / (1 + exp(-2 * KC * A)) + prBias
%   t = A / KC * tanh(KC * A) + t_nd
%   where
%       KC = k * (coh + cohBias) + uBias,
%       k = fitted parameter,
%       A = flat bound height,
%       coh = signed coherence,
%       cohBias = motion strength bias or coherence bias,
%       uBias = drift force bias,
%       t_nd = non-decision time,
%       prBias = vertical probability bias.
%
% To apply Flat Bound fit, first create Flat Bound fit options,

fitOptions = flatBoundFitOptions

% This options contains all the necessary fields to supply raw fitting
% data, fitting and plot options. As shown, 'fitOptions' contains a field 
% 'thetaKey' to identify the key type of each theta value, which is 
%
% 'kappa','A','tndr','tndl','uBias','cohBias','prBias'.
%
% To fit the data using both rightward- and leftward- non-decision time, 
% the 'theta' field can be set as

fitOptions.theta = [0.457*sqrt(1E3),25/sqrt(1E3), 0.3,0.3, NaN,NaN,NaN]

% Please note that any theta parameter that is set to NaN, means that this
% theta parameter is not fitted. In the above case, only 'kappa', 'A',
% 'tndr' and 'tndl' are fitted. 'kappa' & 'A' values are required for all
% fitting.
%
% By default, Flat Bound fit will exclude data if the number of choice is
% less than 'minorRTCriteria' since the variance is not accountable. The
% default 'minorRTCriteria' is 10. If user choose not to reject minor RT,
% 'isRejectMinorRT' field can be set to false.
%
% For this simplest case, Flat Bound fit can be performed by,

[thetaFit,err,exitflag,output,fitOptions] = flatBoundFit(data,fitOptions)

% This call will condense the raw data first, calcute the reaction time for
% each signed coherence, then fit using default FMINSEARCH function and
% outputs
%
%   thetaFit = fitted theta value vector,
%   err = the value of the objective function of fminsearch,
%   exitflat = the exit condition of fminsearch,
%   output = standard output structure of fminsearch,
%   fitOptions = updated fit options with condensed data.
%
% By default, fitting will be performed on condensed data without error RT.
% To fit data with error RT, set 'isFitErrorRT' & 'isPlotErrorRT' field to 
% true. User can also supply fminsearch options by setting fminsearchOptions
% field. To plot the fitting result, call the designed Flat Bound fit 
% function,

fh3 = flatBoundFitPlot(thetaFit,fitOptions)

% Both FLATBOUNDFIT & FLATBOUNDFitPLOT function support the typical Matlab 
% way to change options field value. For example, to plot without the error 
% bar, the plot function can be called as,

fh4 = flatBoundFitPlot(thetaFit,fitOptions,'isPlotErrorBar',false)

% To fit using combined reaction time, simply set the fitting parameter 
% 'theta' field to
%
%   [0.457*sqrt(1E3),25/sqrt(1E3), 0.3,NaN, NaN,NaN,NaN].
%
% and then call the fitting again. 'isFitErrorRT' & 'isPlotErrorRT' fields
% are only available to fitting using both righward- & leftward- reaction
% times. 'isFitCombinedRT' & 'isPlotCombinedRT' fields are used by internal 
% functions to indicate the fitting type. User should not set them.
%
% If fitting using Logit derived coherence bias, set 'ldcb'. For more
% options to fit and explore the Flat Bound fit function, please see also
% FLATBOUNDFIT, FLATBOUNDFITOPTIONS and FLATBOUNDFITPLOT.

%% Monte Carlo Simulation.
% As described in book "Bayesian Brain: Probabilistic Approaches to Neural 
% Coding" edited by Kenji Doya, Shin Ishii, Alexandre Pouget and Rajesh PN.,
% Rao, Chapter 10, "The speed and accuracy of a simple perceptual decision:
% A mathematical primer", the Diffusion-To-Bound framework can be used to
% simulate the accumulation of evidence and the process of decision making. 
% In a similar way, Monte Carlo Method can be used to do this work.
% Furthermore, Monte Carlo Method can be utilized to generate 'fake' Choice
% & Reaction Time data, which else could be fitted by Logit, Weibull & Flat
% Bound Fit to verify the correctness of this method.
%
% The drift term (mu) is
%   kappa*(scoh + cohBias) + uBias
% and the standard deviation term is
%   sqrt(sigma^2 + bSigma * abs(scoh))
%       cohBias = coherence bias or motion strength bias,
%       uBias = drift force bias,
%       sigma & bSigma are 2 parameters to estimate the standard deviation,
%       typically choose 1.0*sqrt(1E3) and 0.5*sqrt(1E3).
% See also MCFIT for more info.
%
% To start with Monte Carlo simulation, first specify the list of signed 
% coherence for generating simulation data as a column vector, 

scoh = [-0.512,-0.256,-0.128,-0.064,-0.032,...
    0,...
    0.032,0.064,0.128,0.256,0.512]'

% Then create Monte Carlo simulation options structure and set 'theta' field
% values. 'theta' field specifies the parameters to calculate the drift term 
% and the standard deviation term. 'thetaKey' field explains the meaning of
% each 'theta' parameter and are,
%
% 'kappa','cohBias','uBias','sigma','bSigma','tndr','tndrsd','tndl','tndlsd'

simOptions = mcSimOptions
simOptions.theta = [0.457*sqrt(1E3), 0.0,0.0, 1.0*sqrt(1E3),0.5*sqrt(1E3), 0.3,0.0, 0.3,0.0]

% The up- & lower- boundary profiles also need to be specified. 6
% predefined boundary profiles are made available and listed in
% GETPROFILEFCN. User can either specify the boundary profile by assigning
% either profile name or function handle. Function handle is in the format
% of @(b,t)(function), where b is a 5-member vector and t is time vector.
% To specify b parameter, set up- & lower- boundary parameter field as
% following. Lower boundary profile is inverted inside MCSIM function to
% create symmetric profile.

simOptions.upBoundaryProfile = 'flat';
simOptions.upBoundaryParameter = [25/sqrt(1E3), 0.0,0.0, 0.0,0.0];

simOptions.lowerBoundaryProfile = 'flat';
simOptions.lowerBoundaryParameter = [25/sqrt(1E3), 0.0,0.0, 0.0,0.0];

% Finally, start Monte Carlo Simulation by

[simData,simOptions] = mcSim(scoh,simOptions)

% In certain situations, signed coherence is only known to be within certain 
% range and distribution. In this case, signed coherence has to be calculated 
% dynamically by supplying a function handle. User is responsible to create 
% this function handle. A sample function handle SCOHDYNSAMPLE was created 
% as example. User also needs to set the total number of simulation trials,

simOptions.trials = 1.0E4

% This 'trials' field now stands for the total number of trials, instead of
% trials per signed coherence used in previous Monte Carlo simulation case.
% To apply Monte Carlo simulation for this dynamic signed coherence, call

[simDataDyn,simOptions] = mcSimDyn(@scohDynSample,simOptions)

% 'simDataDyn' is the calculated simulation data. As shown, signed
% coherence is now dynamic value and no longer distinct repeated value.

%% Diffusion-to-Bound & Monte Carlo fit.
% As mentioned above, both Diffusion to Bound method and Monte Carlo Method
% can be used to fit data. Same as Monte Carlo Simulaiton, the drift term (mu) 
% is
%   kappa*(scoh + cohBias) + uBias
% and the standard deviation term is
%   sqrt(sigma^2 + bSigma * abs(scoh))
%       scoh is signed coherence,
%       cohBias is coherence bias or motion strength bias,
%       uBias is drift force bias,
%       sigma & bSigma are 2 parameters to estimate the standard deviation,
%           typically choose 1.0 and 0.5.
% To set up Diffusion to Bound fit, first create fitting options variable,

fitOptions = dtbMCOptions

% This 'fitOptions' variable is used by both Diffusion to Bound and Monte
% Carlo simulations, and contains several seperate groups of fields.
% The 1st group of fields, 'thetaKey', 'thetaUpLimit', 'theta' &
% 'thetaLowerLimit', specifies key of each parameter, upper & lower limit,
% and initial parameter value to calculate drift term and standard deviation 
% term. DTBMCOPTIONS annotates these fields in detail and the keys of these 
% theta parameters are repeated again,
%
% 'kappa',
% 'cohBias','uBias', 'sigma','bSigma',
% 'tndr','tndrsd', 'tndl','tndlsd',
% 'y0'
%
% These fields can be assigned as following. If the upper limit equals the
% corresponding lower limit and initial fitting value, the optimization
% method will treat this parameter value as fixed and will not optimize it,
% which is the way used here to avoid fitting known parameter. 

fitOptions.thetaUpLimit = [0.8149*sqrt(1E3),...
    0.0,0.0, 1.0*sqrt(1E3),0.0,...
    3.6170,0.5000, 3.6170,0.5000,...
    0.0]
fitOptions.theta = [0.2716*sqrt(1E3),...
    0.0,0.0, 1.0*sqrt(1E3),0.0,...
    0.2794,0.1190, 0.2794,0.1190,...    
    0.0]
fitOptions.thetaLowerLimit = [0.0,...
    0.0,0.0, 1.0*sqrt(1E3),0.0,...
    0.0,0.01, 0.0,0.01,...
    0.0]

% The 2nd group of fields, 'upBoundaryProfile', 'upBoundaryUpLimit',
% 'upBoundaryParameter' & 'upBoundaryLowerLimit', describes the details of
% upper fitting boundary. 'upBoundaryProfile' can be a predefined profile
% string name or function handle in GETFPROFILEFCN, or any user supplied
% function handle in the same format as samples in GETPROFILEFCN.
% 'upBoundaryUpLimit', 'upBoundaryParameter' & 'upBoundaryLowerLimit'
% specifies the up limit, initial fitting value and lower limit of 'b' parameter
% vector of function handle. There is also 3 corresponding fields for the
% lower boundary profile.

fitOptions.upBoundaryProfile = 'flat'
fitOptions.upBoundaryUpLimit = [Inf, 0.0,0.0,0.0,0.0]
fitOptions.upBoundaryParameter = [45.7993/sqrt(1E3), 0.0,0.0,0.0,0.0]
fitOptions.upBoundaryLowerLimit = [0.0, 0.0,0.0,0.0,0.0]

fitOptions.lowerBoundaryProfile = 'flat'
fitOptions.lowerBoundaryUpLimit = [Inf, 0.0,0.0,0.0,0.0]
fitOptions.lowerBoundaryParameter = [45.7993/sqrt(1E3), 0.0,0.0,0.0,0.0]
fitOptions.lowerBoundaryLowerLimit = [0.0, 0.0,0.0,0.0,0.0]

% 'fitType' field is used by internal functions to specify whether DTB
% or Monte Carlo method is chosen for fitting. No need to set
% this field. 'optMethodList' field lists all the available optimization
% methods. User can set 'optMethod' field to one of them for their own favor. 
% By default, if not setting, fminsearchbnd optimization method will be
% used for optimization. 'fmincon' is the recommended optimization method. 
% User can start with 'fmincon' and then try other optimization methods.

fitOptions.optMethod = 'fmincon'

% 'dt' field sets the time step size in fitting. By default, the time step
% size is 1.0E-3 second. It is better to use a time step not coarser that 
% this value. 'tMax' is the maximum time that fitting does. 'rngSeed' is the
% random number generator seed for repeatable random number generation. 
% 'NumWorkers' is the number of processors that parallel computation utilizes. 
% It is recommended to use neither more than the maximum processors that 
% each optimization method utilizes, nor more than the total processors the 
% computer has. If not set, it will use the maximum available processors.

fitOptions.dt = 1.0E-3
fitOptions.tMax = 5.0

% By default, variable duration choice fitting is not enabled. To do 
% variable duration choice fitting, set

fitOptions.isChoiceVariableDuration = true;

% If preferring to plot the iteration results during iteration, 'isPlotIter' 
% field can be set to true. To start DTB fit, call DTBFIT function direcly 
% by supplying fitting data and fitting options. The fitting result is 
% returned as a structure containing fitted parameter and standard returns
% by optimization methods.

fitResult = dtbFit(data,fitOptions)

% When using FP4 (Chang-Cooper Method to solve Fokker-Planck equation) to
% fit data, time step needs to be smaller. Call FPFIT function as
% following.

fitOptions.dt = 0.5E-3
fitResult = fpFit(data,fitOptions)

% For Monte Carlo fit, simply set up the common fit options fields as above.
% Or just reuse the above fitting options. Notice that, 'y0' parameter in 
% 'theta' field of fitOptions is only valid for DTB fit, not valid for 
% Monte Carlo fit. When performing Monte Carlo fitting, simply ignore this
% 'y0' parameter and set its value to 0. Let's set the optimization method
% to 'fmincon',

fitOptions.optMethod = 'fmincon'

% Then start Monte Carlo fit and call MCFIT function,

fitResult = mcFit(data,fitOptions)

% Monte Carlo method can be efficiently parallerized by Nvidia GPU accelerator 
% supported by Matlab. GPU memory is very limited and can be quickly
% exhausted. It is good to start with coarser time step and shorter
% durance. And then gradually use finer time stop and longer durance.

fitOptions.dt = 5.0E-3
fitOptions.tMax = 3.0

% Since GPU memory is small, only 1 computing process is allowed by default
% to avoid draining out GPU memory. Set 'isUseGPU' field to true.

fitOptions.NumWorkers = 1
fitOptions.isUseGPU = true

% Finally, start Monte Carlo fitting with GPU,

fitResult = mcFit(data,fitOptions)







