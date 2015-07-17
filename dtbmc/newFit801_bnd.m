function [logl,uscoh,t,P,pcorr] = newFit801_bnd(theta,data,opt)
%NEWFIT801_BND main function to perform DTB & Monte Carlo fit
%
%   This function returns the -logL value of signed coherence data for either 
%   DTB or Monte Carlo fit.
%
%   [logl,uscoh,t,P,pcorr] = newFit801_bnd(theta,data,opt)
%
%   where      
%       theta is a list of fitting parameters,
%       data is the input data for fit which is in the form of [scoh,choice,rt],
%       and,
%       logl is -1.0 * log-likelihood,
%       uscoh is unique signed coherence level,
%       t is time axis data,
%       P is data structure for dtb_analytic,
%       pcorr is ?.

% To clean up:
% The model assumes that all information is used to make the choice (even
% though a perturbation in the reflex experiment occurs earlier). Note that
% there is no RT in the data. The main return argument, err, is to be
% minimized for fitting. The other arguments are
% x  vector of unique values in column 1 of data (tyPcally signed coh)
% t  time vector (note units are sec)
%    The other arguments are indexed wrt x and distributions are length(t)

%   version 1.0
%   2014 revised basing on Shin Kira version, added new features, renamed file to 
%   newFit801_bnd.m

% Copyright Shadlen Lab 2015.


% Load 'theta' group of fitting parameters.
g1 = theta(1:10);
kappa = g1(1);
cohBias = g1(2);
uBias = g1(3);
sigma = g1(4);
bSigma = g1(5);
tndr = g1(6);
tndrsd = g1(7);
tndl = g1(8);
tndlsd = g1(9);
y0 = g1(10);

scoh = data(:,1); % [-1,1]
ci  = data(:,2); % 0 for leftward choice and 1 for rightward choice.
rt  = data(:,3); % Unit of second.
choice = logical(ci);
uscoh = unique(scoh)'; % Unique signed coherence.

if isempty(opt.dt)
    dt = 0.5E-3; % Unit of second.
else
    dt = opt.dt;
end

if isempty(opt.tMax)
    tmax = 5; % Unit of second.
    % tmax = max(rt)+0.3;
else
    tmax = opt.tMax;
end

t = (0:dt:tmax)';

% load the 2nd group of fitting parameters, calculate up-boundary profile
b = theta(11:15);
Bup = feval(opt.upBoundaryProfile,b,t');
% Bounds stop collapsing when the bounds become less than 0.1% of initial height
if isnan(b(1))
    error('newFit801_bnd:bup(1) must be provided.');
end

Bup(Bup <= b(1)*1e-3, 1) = b(1) * 1e-3;

% load the 2nd group of fitting parameters, calculate lower-boundary profile
b = theta(16:20);

if isempty(b) % symmetrical boundary
    Blow = -Bup; 
else % non-symmetrical boundary
    Blow = feval(opt.lowBoundaryProfile,b,t');
    if isnan(b(1))
        error('newFit801_bnd:blo(1) must be provided.');
    end
    
    Blow(Blow <= b(1)*1e-3, 1) = b(1) * 1e-3;
    Blow = -1.0 * Blow; % inverse, make defining boundary profile easier
end

uv = kappa*(uscoh + cohBias) + uBias; % Drift term (mu)
dfu = sqrt(sigma^2 + bSigma * abs(uscoh)); % diffusion term without normrnd()
rngSeed = opt.rngSeed;
useGPU = opt.isUseGPU;

%TODO: output P must be wrappered into the exact same format.
switch opt.fitType
    case 'Monte Carlo'
        P = mckernel(uv,t,Bup,Blow,dfu,rngSeed,useGPU); % Monte Carlo Fit
        
    case 'DTB'
        if 0 %all(Bup==Bup(1)) && all(Blo==Blo(1)) %we can go analytic
            % P = analytic_dtb(uv,t,Bup(1),Blo(1),y0);
        else  % or need to sikappalate with fft method
            % md=max(abs(uv));
            % sm=(md*dt+sqrt(dt)*4); %DW to test sense of this
            % y=linspace(min(Blo)-sm,max(Bup)+sm,512)';
            nBin = 2^9;
            y = determineY(kappa, uscoh, dt, Bup, Blow, nBin);
            yinit = 0*y;
            
            i1 = find(y>=y0, 1,'first');
            i2 = find(y<=y0, 1,'last');
            if i1 == i2
                yinit(i1)=1;
            else
                w2=abs(y(i1)-y0);
                w1=abs(y(i2)-y0);
                
                w1=w1/(w1+w2);
                w2=(1-w1);
                yinit(i1)=w1;
                yinit(i2)=w2;
            end
        
            % P = spectral_dtb(uv,t,Bup,Blo,y,yinit); % original version
            
            notabs_flag = opt.isChoiceVariableDuration;
            useGPU = false;
            
            % antialiased version spectral_dtb
            if any(abs(dfu) > eps)
                P = spectral_dtbAA(uv,t,Bup,Blow,y,yinit,notabs_flag,useGPU,dfu,rngSeed);
            else
                P = spectral_dtbAA(uv,t,Bup,Blow,y,yinit,notabs_flag,useGPU);
            end
        end
        
    case 'FP4'
        
        
    otherwise
        error('No predefined fit type matched.');
        
end


if ~opt.isChoiceVariableDuration    
    for i = 1:length(scoh)
        Ic = uscoh == scoh(i); % index for the coherence
        It = find(t>=rt(i), 1); % index for the time (end of viewing time or rt)
        
        % distribution of tnd going from rt to zero
        if scoh(i) >= -1.0*cohBias
            tnd = tndr;
            tnd_sd = tndrsd;
        else
            tnd = tndl;
            tnd_sd = tndlsd;
        end
        
        r = normpdf(rt(i)-(1:It)*dt,tnd,tnd_sd)*dt;
        
        % reaction time
        p_up(i) = P.up.pdf_t(Ic,1:It) * r'; %#ok<AGROW>
        p_lo(i) = P.lo.pdf_t(Ic,1:It) * r'; %#ok<AGROW>
    end
    
    % ensure that probabilites lie between eps and 1-eps
    p_up = clip(p_up,eps,1-eps);
    p_lo = clip(p_lo,eps,1-eps);
    pPred = p_up.*choice'+ p_lo.*~choice';
    pcorr = p_up'.*(scoh>0)+ p_lo'.*(scoh<0);
    logl = -sum(log(pPred));

else
    % Choice Variable Duration fit error. Error calculation method ported
    % from choiceVdurFromDiffusion_err.m.
    
    % Variable name 'pPred' is chosen for consistence and is equavalent to
    % variable 'p' in choiceVdurFromDiffusion_err.m.
    pPred = zeros(size(scoh,1),1);
    
    for i = 1:length(scoh)
        Ic = find(uscoh == scoh(i)); % index for the coherence
        It = find(t>=rt(i), 1); % index for the time (end of viewing time or rt)
                        
        if choice(i) == 1
            pPred(i) = P.notabs.pos_t(Ic,It) + P.up.cdf_t(Ic,It);
        else
            pPred(i) = 1 - (P.notabs.pos_t(Ic,It) + P.up.cdf_t(Ic,It));
        end
    end
        
    pcorr = NaN;
    logl = -sum(log(pPred));
end

fprintf('err= %.3f\n',logl);


