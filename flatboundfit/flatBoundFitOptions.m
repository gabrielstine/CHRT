classdef flatBoundFitOptions
    %FLATBOUNDFITOPTIONS
    
    %   Copyright Jian Wang 2015
    
    
    properties (GetAccess = public, SetAccess = public)
        theta = NaN(1,7);
        thetaKey = {'kappa',...
            'A',... % Boundary height
            'tndr',... % Rightward non-decision time or Combined non-decision time
            'tndl',... % Leftward non-decision time
            'uBias',... % Drift Force bias
            'cohBias',... % Motion Strength Bias
            'prBias' % Vertical Probability Bias
            };
        thetaFit;
        
        condensedData;
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

