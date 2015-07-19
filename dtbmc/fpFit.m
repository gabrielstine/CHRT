function fitResult = fpFit(data,fitOptions,varargin)
%FPFIT apply Fokker Planck finite differnece simulation to fit 1D signed 
%coherence data
%
%   fitResult = fpFit(icd,varargin)
%
%   This fit function is a wrapper of DTBMC and uses FP4 for
%   Fokker Planck finite difference simulation.
%
%   See also DTBMC, DTBFIT, MCFIT, DTBMCOPTIONS.

%   Copyright 2014 Jian Wang

fitOptions.fitType = 'FP4';
% DTBMC kernel simulates DTB, MC & FP4.
fitResult = dtbmc(data,fitOptions,varargin{:});
