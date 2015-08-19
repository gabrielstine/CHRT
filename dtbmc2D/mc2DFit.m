function fitResult = mc2DFit(data,fitOptions,varargin)
%MC2DFIT apply Monte Carlo method to fit 1D signed coherence data
%
%   fitResult = mc2DFit(data,fitOptions,varargin)
%
%   This fit function is a wrapper of DTBMC and uses MCKERNEL2D for Monte Carlo
%   Simulation.
%
%   See also DTBMC, DTBMC2DOPTIONS.

%   Copyright 2015 Jian Wang

fitOptions.fitType = 'Monte Carlo 2D';
% DTBMC kernel simulates both DTB & MC.
fitResult = dtbmc2D(data,fitOptions,varargin{:});
clear mckernel2D;

