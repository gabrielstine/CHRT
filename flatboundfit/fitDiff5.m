function err = fitDiff5(theta,data,thetaFlag)
%FITDIFF5 main function to calculate Flat Bound fitting error
%   err = fitDiff5(theta,data,thetaFlag) calculate Flat Bound fitting error
%   where 
%       theta is fitting parameter vector,
%       data is fitting data vector,
%       thetFlag is a bool vector that shows which 'theta' parameters are 
%       provided for fitting,
%       and
%       err returns the calculated fitting error.
%
%   See also CALDIFF5.

%   Copyright Shadlen Lab 2015.

cohs   = data(:,1);
t1_obs = data(:,2);
t1_se  = data(:,3);
t2_obs = data(:,4);
t2_se  = data(:,5);
n1_obs = data(:,6);
n_total = data(:,7);

[t1_pred,t2_pred,p_pred] = calcDiff5(cohs,theta,thetaFlag);

% Cost function of RT.
N = size(data,1);

t1_ind = ~isnan(t1_pred) & ~isnan(t1_se);
err = sum((t1_obs(t1_ind) - t1_pred(t1_ind)).^2./(2.*t1_se(t1_ind).^2)) ...
    + N./2.*log(2.*pi) + sum(log(t1_se(t1_ind))); % TODO: remove t1_se err term

t2_ind = ~isnan(t2_pred) & ~isnan(t2_se);
err = err + sum((t2_obs(t2_ind) - t2_pred(t2_ind)).^2./(2.*t2_se(t2_ind).^2)) ...
    + N./2.*log(2.*pi) + sum(log(t2_se(t2_ind))); % TODO: remove t2_se err term 
    
% Cost function of choice probability.
err = err - sum(gammaln(n_total+1) - gammaln(n1_obs+1) - gammaln(n_total-n1_obs+1)...
    + n1_obs.*log(p_pred+eps) + (n_total-n1_obs).*log(1-p_pred+eps));




