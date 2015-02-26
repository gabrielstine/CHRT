function fitOptions = initParallel(fitOptions)
%INITPARALLEL initialize parallel computing environment

% Copyright 2014 Jian Wang

if fitOptions.isUseGPU
    if gpuDeviceCount < 1
        error('No GPU accelerator detected.');
    end
    
    fitOptions.NumWorkers = 1; % Use 1 process by force.
    
    % To do: handle GPU hardware computing capability.
end


ph = gcp('nocreate'); % Get current pool handle.

if isempty(ph)
    currentWorkers = 0;
else
    currentWorkers = ph.NumWorkers;
end

if isempty(fitOptions.NumWorkers)
    thisCluster = parcluster();
    newWorkers = thisCluster.NumWorkers; % NumWorkers of default local profile.
else
    newWorkers = fitOptions.NumWorkers;
end

if currentWorkers ~= newWorkers
    delete(ph);
    ph = parpool(newWorkers); % Create new pool using local profile.
    
    % To do: handle cluster profile.
end


