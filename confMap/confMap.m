function [pSet,M] = confMap(c,w,varargin)
% constructs a confidence map based on prior over c whose weights are in w,
% the prior over c. 
% see theory/test_confMap.m



% history
% 3/4/15 mns initiated after finding a bug in confidenceMap
% 6/6/2015 mns added kappa

%% Parse input
pSet = inputParser;


 
% addParamValue(pSet,'t',@iscolumn);
addParamValue(pSet,'dt',0.001);
addParamValue(pSet,'tmin',0);  % earliest time on graph
addParamValue(pSet,'tmax',3);  % earliest time on graph
addParamValue(pSet,'ymax',1,@(x)x>0);
addParamValue(pSet,'ny',100); % steps in the ygrid
addParamValue(pSet,'kappa',1);

% addParamValue(pSet,'showGraph',true,@islogical);
% addParamValue(pSet,'newFig',false,@islogical);
% addParamValue(pSet,'title','');

parse(pSet,varargin{:});

%% pull out some variables
if isscalar(c)
    c = [-c;c];
    w = .5*[1;1];
    warning('applying symmetric prior for c')
end

if ~(iscolumn(c) && iscolumn(w)) | ~all(size(c)==size(w))
    error('c and w must be same size column vectors');
end
if sum(w)~=1
    w = w./sum(w);
end

if ~all(ismember([1,-1],unique(sign(c))))
    % user supplied pos or neg coh ±0
    c0=c;
    [c,I] = sort([c;-c(c~=0)]);
    w = [w;w(c0~=0)];
    w = w(I)./sum(w);
    warning('applying symmetric prior for c')
end

t = (pSet.Results.tmin:pSet.Results.dt:pSet.Results.tmax)';
y = linspace(0,pSet.Results.ymax, pSet.Results.ny + 1)';
y = sort([y;-y(y>0)]);


%% compute the map
% make the grid
[T,Y] = meshgrid(t,y);
PC = zeros(size(Y)); % initialize correct
PE = zeros(size(Y)); % initialize error

for i = 1:length(c)
    p = normpdf(Y,kappa*c(i)*T,sqrt(T));
    if c(i)>0
        PC(Y>0) = PC(Y>0) + w(i)*p(Y>0);
        PC(Y==0) = PC(Y==0) + w(i)*0.5;
        PE(Y==0) = PE(Y==0) + w(i)*0.5;
        PE(Y<0) = PE(Y<0) + w(i)*p(Y<0);
    elseif c(i)==0
        PC = PC + w(i)*0.5;
        PE = PE + w(i)*0.5;
    else
        PE(Y>0) = PE(Y>0) + w(i)*p(Y>0);
        PE(Y==0) = PE(Y==0) + w(i)*0.5;
        PC(Y==0) = PC(Y==0) + w(i)*0.5;
        PC(Y<0) = PC(Y<0) + w(i)*p(Y<0);
    end
end



%% 
        
    

%% return
M.c = c;
M.w = w;
M.k = kappa;
M.t = t;
M.y = y;
M.PC = PC;
M.PE = PE;
M.CMAP = log(PC./PE);









