function [F f] = rayleighCDF(t,dt,fshift,sigma)
% Reyleigh distribution for time-out.

t = t*10^-3;
f = (t-fshift).*exp(-(t-fshift).^2./(2.*sigma.^2))./(sigma.^2);
f(f(:)<0) = 0;
F = cumsum(f).*dt.*10^-3;