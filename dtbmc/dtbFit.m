function fitResult = dtbFit(data,fitOptions,varargin)
%DTBFIT apply Diffusion to Bound simulation to fit 1D signed coherence data
%
%   fitResult = dtbfit(icd,varargin)
%
%   This fit function is a wrapper of DTBMC and uses SPECTRAL_DTBAA for
%   Diffusion to Bound simulation.
%
%   See also DTBMC, MCFIT, DTBMCOPTIONS.

%   Copyright 2014 Jian Wang

fitOptions.fitType = 'DTB';
% DTBMC kernel simulates both DTB & MC.
fitResult = dtbmc(data,fitOptions,varargin{:});
