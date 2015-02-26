function fh = weibullFitPlot(alpha,beta,pcData) 
%WEIBULLFITPLOT plot Weibull Fit data
% fh = weibullFitPlot(alpha,beta, pcData)
%   where
%       alpha is threshold parameter, 
%       beta is slope parameter,
%       pcData is the proportional correct choice data,
%       and
%       fh is the figure handle.

%   Copyright Jian Wang 2015

% Calculate weibull fit data.
fWeib = @(q,x) 1 - 0.5 * exp(-(x/q(1)).^q(2));
xrange = 0.55;
xf = linspace(-0.1*xrange,1.0*xrange,101);
yf = real(fWeib([alpha,beta],xf));

% Plot ucoh vs. pccho.
fh = figure;
hold on;
plot(pcData(:,1),pcData(:,2),'LineStyle','none','Marker','o',...
    'MarkerFaceColor','b');
plot(xf,yf,'m');
hold off;

set(gca,'xlim',xrange*[-0.1,1.0]);
title('Weibull Fit','FontSize',14,'FontWeight','bold');
xlabel('Motion strength');
ylabel('Proportion correct choices');

end
