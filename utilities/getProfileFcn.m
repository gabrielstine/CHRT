function fh = getProfileFcn(profname)
%GETPROFILEFCN return predefined boundary profile handle
%   fh = getProfileFcn(profname)
%       profname can be profile name string or a function handle.
%       fh is the profile function handle.

%   Copyright 2015 Jian Wang

if ischar(profname)
    profname = lower(profname);
    switch profname
        case 'linear'
            fh = @(b,t) cat(2,repmat(b(1),1,length(t(t<=b(3)))),...
                b(1)-b(2)*(t(t>b(3))-b(3)));
        case 'quadratic'
            fh = @(b,t) cat(2,repmat(b(1),1,length(t(t<=b(3)))),...
                b(1)-b(2)*(t(t>b(3))-b(3)).^2);
        case 'exponential'
            fh = @(b,t) cat(2,repmat(b(1),1,length(t(t<=b(3)))),...
                b(1)*exp(-b(2)*(t(t>b(3))-b(3))));
        case 'logit'
            fh = @(b,t) b(1)./(1+exp(b(2)*(t-b(3))));
        case 'hyperbolic'
            fh = @(b,t) cat(2,repmat(b(1),1,length(t(t<=b(3)))),...
                b(1)./(1+b(2)*(t(t>b(3))-b(3))));
        otherwise  % 'flat'
            fh = @(b,t) repmat(b(1),1,length(t));
    end
elseif isa(profname,'function_handle')
    
else
    error('Profile name needs to be either a string or a function handle.');
end
