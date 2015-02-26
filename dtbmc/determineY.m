function y = determineY(k, coh, dt, ub, lb, ny, sig, marginFactor) % , granularity)
% y = determineY(k, coh, dt, ub, lb, ny, sig, marginFactor) % , granularity)
%
% Determines y for spectral_dtb type functions.
%
% Required arguments:
%
% k  : multiplies coh
% coh: vector of coherences
% dt : size of the time step
% ub : upper bound. Can be a vector.
% lb : lower bound. Can be a vector.
%
% Optioinal arguments:
%
% ny : minimum length of y, say, 2^12. 
% sig: standard deviation.
% marginFactor: minimum margin as a factor of sig.
% granularity: minimum granularity that divides sig. % unused.

% Copyright Shadlen Lab 2015.

if ~exist('ny' , 'var') || isempty(ny), ny  = 2^12; end
if ~exist('sig', 'var') || isempty(sig), sig = 1; end
if ~exist('marginFactor', 'var') || isempty(marginFactor), marginFactor = 4; end
% if ~exist('granularity', 'var')  || isempty(granularity), granularity = 20; end

eCoh = k*coh*dt;
eSig = sig*sqrt(dt);

% minY = min(min(lb), min(eCoh)) - marginFactor*eSig;
maxY = max(max(ub), max(eCoh)) + marginFactor*eSig;
minY = -maxY;

dy   = min((maxY-minY)/ny); % , eSig/granularity);
ny   = pow2(nextpow2((maxY-minY)/dy));

y    = linspace(minY,maxY,ny)';
