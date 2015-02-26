classdef mcSimOptions
    %MCSIMOPTIONS Monte Carlo Simulation Options   
    
    %   Copyright 2014 Jian Wang    
    
    properties (GetAccess = public, SetAccess = public)                
        dt@double = 0.1E-3; % Time step size (unit of second).
        tMax@double = 5.0; % Simulation time length.
        
        rngSeed; % Random number generator seed.
                
        % For simulation case by providing unique signed coherence,
        % this 'trials' number means the number of trials per signed coherence.
        % For simulation case by providing signed coherence function,
        % the signed coherence is calculated out from a nomral distribution
        % function and this 'trials' number means the total number of
        % simulations. Each calculated signed coherence is simulated single 
        % once.
        trials@double = 100;
                
        theta = zeros(1,9); % Default values are zeros.
        % sigma & bSigma are 2 parameters to calculate the standard deviation 
        % of drift term as in sqrt(sigma^2 + bSigma * abs(scoh)). See also
        % MCSIM.
        thetaKey = {'kappa',... 
            'cohBias',... % Coherence bias or Motion Strength Bias.
            'uBias',... % Drift Force bias.
            'sigma',...
            'bSigma',...
            'tndr',... % Rightward non-decision time.
            'tndrsd',... % Rightward non-decision time standard deviation.
            'tndl',... % Leftward non-decision time.
            'tndlsd' % Leftward non-decision time standard deviation.
            };
        
        upBoundaryProfile; % Up-bounday profile name or function handle.
        upBoundaryParameter = zeros(1,5); % Default values are zeros.
        
        lowerBoundaryProfile; % Lower-bounday profile name or function handle.
        lowerBoundaryParameter = zeros(1,5); % Default values are zeros.                                        
    end
    
end