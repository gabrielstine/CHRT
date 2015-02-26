classdef dtbMCOptions
    %DTBMCOPTIONS Diffusiont-To-Bound & Monte Carlo fitting options.
    %
    %   thetaKey - Cell row vector of strings to identify each theta fitting    
    %       parameter.
    %   thetaUpLimit - Row vertor of upper bound of theta fitting
    %       parameters.        
    %   theta - Row vector of initial theta fitting parameters.
    %   thetaLowerLimit - Row vector of lower bound of theta fitting
    %       parameters.
    %   
    %   upBoundaryProfile - String keyword (flat,linear,quadratic,exponential,
    %       logit,hyperbolic) or function handle to pass the upper boundary 
    %       profile.    
    %   upBoundaryUpLimit - Row vector of upper limit of upper boundary profile 
    %       parameters.
    %   upBoundaryParameter - Row vector of upper boundary profile parameters.
    %   upBoundaryLowerLimit - Row vector of lower limit of upper boundary 
    %       profile parameters.
    %    
    %   lowerBoundaryProfile - String keyword or function handle to pass the
    %           lower boundary profile. The lower boundary profile can be passed 
    %           in the same format as the up boundary profile since it is 
    %           automatically inverted inside.
    %   lowerBoundaryUpLimit - Row vector of upper limit of lower boundary 
    %       profile parameters.
    %   lowerBoundaryParameter - Row vector of lower boundary profile parameters.
    %   lowerBoundaryLowerLimit - Row vector of lower limit of lower boundary 
    %       profile parameters.
    %    
    %   fitType - Specify either DTB or Monte Carlo Fit, used internally by
    %       DTBFIT & MCFIT function.
    %
    %   optMethod - String to specify which optimization method used for 
    %       fit.            
    %   optMethodList - String cells to list all the available optimization
    %       methods.
    %            
    %   dt - Simulation time step size.
    %   tMax - Max simulation time.
    %
    %   rngSeed - Random number generator seed.
    %   NumWorkers - Number of parallel computation threads.
    %   isUseGPU - Whether to use GPU for computing or not.
    %
    %   isPlotIter - Whether plot the iteration result or not.             
                
    %   Copyright 2014 Jian Wang

    properties (GetAccess = public, SetAccess = public)                
                
        thetaKey = {'kappa',...
            'cohBias',... % Coherence bias or motion strength bias.
            'uBias',... % Drift force bias.
            'sigma',... % Parameters to estimate the standard deviation.
            'bSigma',...            
            'tndr',... % Rightward non-decision time.
            'tndrsd',... % Standard deviation of rightward non-decision time.
            'tndl',... % Leftward non-decision time.
            'tndlsd',... % Standard deviation of leftward non-decision time.
            'y0'}; % ???       
        thetaUpLimit@double; % Up limit of theta field.
        theta@double;
        thetaLowerLimit@double; % Lower limit of theta field.
                               
        upBoundaryProfile; % Upbounday profile name or function handle.        
        upBoundaryUpLimit@double; % Up limit of upboundary profile parameter vector.
        upBoundaryParameter@double; % Upboundary profile parameter vector.
        upBoundaryLowerLimit@double; % Lower limit of upboundary profile parameter vector.
        
        lowerBoundaryProfile; % Lower-bounday profile name or function handle.        
        lowerBoundaryUpLimit@double; % Up limit of lower boundary profile parameter vector.
        lowerBoundaryParameter@double; % Lower boundary profile parameter vector.
        lowerBoundaryLowerLimit@double; % Lower limit of lower boundary profile parameter vector.
        
        fitType; % DTB or Monte Carlo Fit.
        
        optMethod; % Optimization method.                
        optMethodList = {'fminsearch',...
            'fminsearchbnd',...
            'fmincon',...
            'globalsearch',...
            'multistart',...
            'patternsearch'}; % List of optimization methods.
        
        dt@double; % Simulation time step size.
        tMax@double % Max simulation time.
        
        rngSeed; % Random number generator seed.        
        NumWorkers@double; % Number of parallel computation workers.                
        isUseGPU@logical = false; % Whether to use GPU for computing or not.
                
        isPlotIter@logical = false; % Plot the iteration result or not.             
    end                                                         
end

