function [thetaFit,err,exitflag,output,fitOptions] = flatBoundFit(data,fitOptions,varargin)
%FLATBOUDFIT apply Flat Bound fit to 1D signed coherence data
%   [thetaFit, fitOptions] = flatBoundFit(data,fitOptions,varargin)
%   where
%       data is 1D signed coherence data,
%       fitOptions is Flat Bound fit options,
%       varargin are optional <field,value> pairs to update fitting options,
%       and
%       thetaFit is fitted theta vector,
%       err is the fitting value of the objective function,
%       exitflat is the exit condition of fminsearch,
%       output contains the information about fminsearch optimization,
%       fitOptions is updated fitOptions containing condensed data.
%
%   See also flatBoundFitOptions, fitDiff5, FMINSEARCH.

%   Copyright Jian Wang 2014

if nargin < 2
    error('Not enough inputs.');
end

if nargin > 2
    fitOptions = updateOptions(fitOptions,varargin{:}); % Update options
end

% Check must-be supplied theta parameters
if isnan(fitOptions.theta(1))
    error('Invalid Kappa value.');
end

if isnan(fitOptions.theta(2))
    error('Invalid A value.');
end

if isnan(fitOptions.theta(3))
    error('Invalid tndr value.');
end

scoh = unique(data(:,1)); % Number of unique signed coherence
condensedData = zeros(length(scoh),length(fitOptions.condensedDataKey)); 
ldcb = -1.0 * fitOptions.ldcb; % Logit derived coherence bias

for i = 1:length(scoh) % Condense data
    s = data(:,1) == scoh(i);
    nrcho = sum(data(s,2) == 1); % Number of rightward choice
    rtr = data((s & (data(:,2)==1)), 3); % Rightward RT
    rtl = data((s & (data(:,2)==0)), 3); % Leftward RT
    rtrm = mean(rtr);
    rtrse = std(rtr) ./ sqrt(length(rtr));
    rtlm = mean(rtl);
    rtlse = std(rtl) ./ sqrt(length(rtl));
    
    % Exclude data that number of choices are less than "MinorRTCriteria"
    % since the variance is not accountable
    if fitOptions.isRejectMinorRT        
        if nrcho < fitOptions.minorRTCriteria
            rtrm = NaN;
            rtrse = NaN;
        end
        
        if sum(s) - nrcho < fitOptions.minorRTCriteria
            rtlm = NaN;
            rtlse = NaN;
        end
    end
    
    % Combine RT to calculate RT_correct
    if scoh(i) > ldcb % Rightward choice is correct
        rtcm = rtrm;
        rtcse = rtrse;
    elseif scoh(i) < ldcb % Leftward choice is correct
        rtcm = rtlm;
        rtcse = rtlse;
    elseif scoh(i) == ldcb % All choice are correct
        rtc = data(s,3);
        rtcm = mean(rtc);
        rtcse = std(rtc) ./ sqrt(length(rtc));
    end
    
    condensedData(i,:) = [scoh(i),rtrm,rtrse,rtlm,rtlse,nrcho,sum(s),rtcm,rtcse];
end

fitOptions.condensedData = condensedData; % Further changes will not be reserved

theta = fitOptions.theta;
if sum(isnan(theta)) == 4
    fitOptions.isFitCombinedRT = true;
    fitOptions.isPlotCombinedRT = true;
elseif sum(isnan(theta)) == 3
    fitOptions.isFitCombinedRT = false;
    fitOptions.isPlotCombinedRT = false;
end

% Prepare data for optimization
if fitOptions.isFitCombinedRT % Case of fitting using only combined RT
    if ~isnan(theta(4))
       error('tndl is provided.'); 
    end
    
    d4f = NaN(length(scoh),7); % [scoh,rtcm,rtcse,NaN,NaN,nrcho,nt]
    d4f(:,1) = condensedData(:,1);
    d4f(:,2:3) = condensedData(:, 8:9);
    d4f(:,6:7) = condensedData(:, 6:7);
else
    if isnan(theta(4))
       error('tndl is not provided.'); 
    end
    
    d4f = condensedData(:,1:7); % use RT with error
    
    if ~fitOptions.isFitErrorRT % use RT without error
        d4f(scoh < ldcb, 2:3) = NaN;
        d4f(scoh > ldcb, 4:5) = NaN;
    end
end

theta = theta(~isnan(theta));

% Optimize using fminsearch
if isempty(fitOptions.fminsearchOptions)
    fminsearchOptions = optimset('fminsearch');
    fminsearchOptions = optimset(fminsearchOptions,'MaxFunEvals',10.^5,'MaxIter',10.^5);
end

[thetaFit,err,exitflag,output] = fminsearch('fitDiff5',theta,fminsearchOptions,d4f,~isnan(fitOptions.theta));
fitOptions.thetaFit = thetaFit;

end