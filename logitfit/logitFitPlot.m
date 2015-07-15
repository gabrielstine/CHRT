function fh = logitFitPlot(prData,b)
%LOGITFITPLOT plot Logit Regression Fit data
%   fh = logitFitPlot(prData,b)
%   where
%       prData is proportion rightward choice data, and
%       b is a vector of logit fit coefficient estimates.

%   Copyright Jian Wang 2015

p = inputParser;
addRequired(p,'prData',@(x) ismatrix(x) && size(x,2) == 2);
addRequired(p,'b',@(x) isvector(x) && length(x) == 2);
parse(p,prData,b);

% Calculate Logit Regression fitting data.
xrange = 0.55;
x = linspace(-1.0*xrange,1.0*xrange,101);
fLogist = @(bf,xf) 1 ./ (1+exp(-bf(1)-bf(2).*xf)); % Logit Fit function
y = fLogist(b,x);

% Plot scoh vs. prcho.
fh = figure;
plot(prData(:,1),prData(:,2),'LineStyle','none','Marker','o','MarkerFaceColor','b');
hold on;
plot(x,y,'m');
hold off;

set(gca,'xlim',xrange*[-1.0,1.0]);
title('Logit Regression Fit','FontSize',14,'FontWeight','bold');
xlabel('Motion strength');
ylabel('Proportion rightward choices');
