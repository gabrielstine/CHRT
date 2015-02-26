function options = updateOptions(options,varargin)
%UPDATEOPTIONS update options <field,value> pair
%   options = updateOptions(options,field1,value1,...) updates options 
%   struct <field,value> pair
%   where
%       options is options class and
%       varargin is <field,value> pairs for update.

%   Copyright Jian Wang 2014

i = 1; % Index of varargin

if rem(nargin-i,2) ~= 0 % Verify <field,value> be pairs
    error('<field,value> must be pairs.');
end

while i < nargin % Loop through <field,value> pairs
    if any(strcmp(varargin{i},fieldnames(options)))
        options.(varargin{i}) = varargin{i+1};
    else
        warning('Invalid option: %s', varargin{i});
    end
    
    i = i+2;
end

