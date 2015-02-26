function fh = flatBoundFitPlot(thetaFit,fitOptions,varargin)
%FLATBOUNDFITPLOT plot Flat Bound fit data
%   fh = flatBoundFitPlot(thetaFit,fitOptions,varargin)
%       thetaFit is the fitted theta parameter vector,
%       fitOptions is Flat Bound fit options variable,
%       varargin is variable-input length argument list to update fitOptions 
%           variable,
%       fh is the flag handle.

%   Copyright Jian Wang 2015.

if nargin > 1
    fitOptions = updateOptions(fitOptions,varargin{:}); % Update options
end

% Calculate Flat Bound fit data.
xrange = 0.55;
xf = linspace(-1.0*xrange,1.0*xrange,101);
[rtrmp,rtlmp,prchop] = calcDiff5(xf,thetaFit,~isnan(fitOptions.theta));

data = fitOptions.condensedData;
% dataKey = obj.dataKey;
% v = values(containers.Map(dataKey,1:length(dataKey)),...
%     {'scoh','rtrm','rtrse','rtlm','rtlse','nrcho','nt','rtcm','rtcse'});

dfp = get(0,'defaultfigureposition');
dfp(4) = dfp(4)*2.2;
fh = figure('Position',dfp); % Create figure with optimal size

% Plot scoh vs. prcho.
subplot(2,1,1);
hold on;
prcho = data(:,6) ./ data(:,7);
plot(data(:,1),prcho,... % Plot experiment data
    'LineStyle','none','Marker','o','MarkerFaceColor','b');
plot(xf,prchop,'m'); % Plot Flat Bound Fit data
hold off;
set(gca,'xlim',xrange*[-1,1]);
title('Flat Bound Fit','FontSize',14,'FontWeight','bold');
xlabel('Motion strength');
ylabel('Proportion rightward choices');

% Plot scoh vs. rt
subplot(2,1,2);
srlg = {}; % List of legends
hold on;
% Plot experiment data
if fitOptions.isPlotCombinedRT % Plot combined RT
    if fitOptions.isPlotErrorBar
        errorbar(data(:,1),data(:,8),data(:,9),...
            'LineStyle','none','Marker','o','MarkerFaceColor','b');
    else
        plot(data(:,1),data(:,8),...
            'LineStyle','none','Marker','o','MarkerFaceColor','b');
    end
    srlg{end+1} = 'RT combined';
else % Plot non-combined RT
    if ~fitOptions.isPlotErrorRT % Screen error RT
        ldcb = -1.0 * fitOptions.ldcb;
        data(data(:,1) < ldcb,2) = NaN;
        data(data(:,1) < ldcb,3) = NaN;
        data(data(:,1) > ldcb,4) = NaN;
        data(data(:,1) > ldcb,5) = NaN;
    end
    
    if fitOptions.isPlotErrorBar
        errorbar(data(:,1),data(:,2),data(:,3),...
            'LineStyle','none','Marker','o','MarkerFaceColor','b');
        errorbar(data(:,1),data(:,4),data(:,5),...
            'LineStyle','none','Marker','o','MarkerFaceColor','m');
    else
        plot(data(:,1),data(:,2),...
            'LineStyle','none','Marker','o','MarkerFaceColor','b');
        plot(data(:,1),data(:,4),...
            'LineStyle','none','Marker','o','MarkerFaceColor','m');
    end
    srlg{end+1} = 'RT right';
    srlg{end+1} = 'RT left';
end

plot(xf,rtrmp,'b',xf,rtlmp,'m'); % Flat Bound fit data
srlg{end+1} = 'Flat Bound fit';
if ~fitOptions.isFitCombinedRT
    srlg{end+1} = 'FlatB fit';
end
hold off;
legend(srlg{:});
set(gca,'xlim',xrange*[-1,1]);
xlabel('Motion strength');
ylabel('Reaction time (second)');
