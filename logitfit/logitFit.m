function [b,dev,stats,prData] = logitFit(data)
%LOGITFIT apply Logit Regression fit to 1D signed coherence data.
%   [b,dev,stats,prData] = logitFit(data)
%   where
%       data is 1D singed coherence data in the form of [signed coherence 
%       (-1.0, 1.0), choice (0/Left, 1/Right), reation time (second)], and
%       b is a vector of coefficient estimates,
%       dev is the deviance of the fit,
%       stats is the structure containing fitting fields,
%       prData is the calculated proportion rightward choice data.
%       
%   See also GLMFIT.

%   Copyright Jian Wang 2014

p = inputParser;
addRequired(p,'data',@(x) ismatrix(x) && size(x,2) == 3);
parse(p,data);
             
scoh = unique(data(:,1)); % Unique signed coherence.
prData = zeros(length(scoh),2);

% Calculate proportion rightward choice rate.
for i = 1:length(scoh)
    s = (data(:,1) == scoh(i));
    prCho = sum(data(s,2) == 1) / sum(s);
    prData(i,:) = [scoh(i),prCho];
end

% Fit condensed signed data using Logit Regression.    
[b,dev,stats] = glmfit(prData(:,1),prData(:,2),'binomial','link','logit');


