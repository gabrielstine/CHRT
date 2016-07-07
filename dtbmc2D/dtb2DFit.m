function fitResult = dtb2DFit(data,fitOptions,varargin)
%DTB2DFIT apply Spectral DTB 2D method to fit 1D signed coherence data
%
%   fitResult = dtb2DFit(data,fitOptions,varargin)
%
%   This fit function is a wrapper of DTBMC and uses SPECTRAL_DTB_2D for 2D
%   fitting.
%
%   See also DTBMC, DTBMC2DOPTIONS.

%   Copyright 2015 Jian Wang

fitOptions.fitType = 'DTB 2D';
% DTBMC kernel simulates both DTB & MC.
fitResult = dtbmc2D(data,fitOptions,varargin{:});


