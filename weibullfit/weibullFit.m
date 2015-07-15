function [alpha,beta,llik,abse,pcData] = weibullFit(data)
%WEIBULLFIT apply weibull fit to 1D unsigned coherence data
%   [alpha,beta,llik,abse,theta,pcData] = weibullFit(data)
%   where
%       data is 1D singed coherence data in the form of [signed coherence 
%       (-1.0, 1.0), choice (0/Left, 1/Right), reation time (second)], and
%       alpha is threshold parameter, 
%       beta is slope parameter,
%       llik is log likelihood of obtaining data given the fit,
%       abse is a 2-vector with standard errors for alpha and beta,
%       pcData is the proportional correct choice data.
%       
%   See also QUICKFIT.

%   Copyright Jian Wang 2014

p = inputParser;
addRequired(p,'data',@(x) ismatrix(x) && size(x,2) == 3);
parse(p,data);
       
% Calculate the proportion of correct choice, RT mean/se using unsigned
% coherence.
ucoh = unique(abs(data(:,1))); % Unique unsigned coherence
pcData = zeros(length(ucoh),3);

for i = 1:length(ucoh)
    if ucoh(i) == 0 % For zero coherence, all choices are correct
        rclc = (data(:,1) == ucoh(i));
        nt = sum(rclc);
        pcCho = 0.5;
    else
        rc = (data(:,1) == ucoh(i)); % Rightward choices
        lc = (data(:,1) == -ucoh(i)); % Leftward choices
        cc = (rc & (data(:,2) == 1)) | (lc & (data(:,2) == 0));
        nt = sum(rc + lc);
        pcCho = sum(cc) / nt;
    end
    pcData(i,:) = [ucoh(i),pcCho,nt];
end

% Fit condensed unsigned data with quickfit.
[alpha,beta,llik,abse,theta] = quickfit(pcData);            

end


