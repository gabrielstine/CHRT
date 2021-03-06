classdef flatBoundFitOptions
    %FLATBOUNDFITOPTIONS
    %   Fitting options used by FLATBOUNDFIT. Read flatBoundFitOptions.m
    %   for more details.
    %
    %   Copyright Jian Wang 2015
    
    
    properties
        theta = NaN(1,7); % Fitting parameters in the order of thetaKey
    end
    
    properties (Constant = true)
        thetaKey = {'kappa',...
            'A',... % Boundary height
            'tndr',... % Rightward non-decision time or Combined non-decision time
            'tndl',... % Leftward non-decision time
            'uBias',... % Drift Force bias
            'cohBias',... % Motion Strength Bias
            'prBias' % Vertical Probability Bias
            };
    end
    
    properties
        thetaFit; % Fitted parameters        
        condensedData;
    end
    
    properties (Constant = true)
        condensedDataKey = {'scoh',...% Unique signed coherence
            'rtrm',... % Mean rightward reaction time
            'rtrse',... % Standard error of mean rightward reaction time
            'rtlm',... % Mean leftward reaction time
            'rtlse',... % Stardard error of mean leftward reaction time
            'nrcho',... % Number of rightward choice
            'nt',... % Number of both choice
            'rtcm',... % Mean correct reaction time
            'rtcse' % Standard error of mean correct reaction time
            };
    end
       
    properties
        ldcb@double = 0; % Logit derived coherence bias.
        
        isRejectMinorRT@logical = true; 
        minorRTCriteria@double = 10;
        
        isFitCombinedRT@logical; % Field used by internal functions to specify
                                 % whether fit using combined RT or not.   
        isPlotCombinedRT@logical;
        
        isFitErrorRT@logical = false;
        isPlotErrorRT@logical = false;
        
        fminsearchOptions;
        
        isPlotErrorBar@logical = true;
    end
    
end

