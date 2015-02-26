function err = psychAnalyticalFit(theta, data)

cohs  = data(:,1);
t_obs = data(:,2);
t_se  = data(:,3);

[t_pred,p_pred] = psychAnalyticalCalc(cohs, theta);

n_obs  = data(:,4);
n = data(:,5);
N = size(data,1);

err = 0;

% Cost function of RT
err = sum((t_obs - t_pred).^2./(2.*t_se.^2)) + N./2.*log(2.*pi) + sum(log(t_se));

% Cost function of choice probability

err = err - sum(gammaln(n+1) - gammaln(n_obs+1) - gammaln(n-n_obs+1)...
    + n_obs.*log(p_pred+eps) + (n-n_obs).*log(1-p_pred+eps));

