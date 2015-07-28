function fitResult = mcFit(data,fitOptions,varargin)
%MCFIT apply Monte Carlo method to fit 1D signed coherence data
%
%   fitResult = mcFit(data,fitOptions,varargin)
%
%   This fit function is a wrapper of DTBMC and uses MCKERNEL for Monte Carlo
%   Simulation.
%
%   See also DTBFIT, DTBMC, DTBMCOPTIONS.

%   Copyright 2014 Jian Wang

fitOptions.fitType = 'Monte Carlo';
% DTBMC kernel simulates both DTB & MC.
fitResult = dtbmc(data,fitOptions,varargin{:});
clear mckernel;

