function D =  mckernel(drift,t,Bup,Blo,y0,sigma,rngseed,notabs_flag,useGPU,trials)
%MCKERNEL kernel of Monte Carlo simulation for MCFIT fitting
%
% D =  mckernel(drift,t,Bup,Blo,y0,sigma,rngseed,notabs_flag,useGPU,trials)
%
% This Monte Carlo Simulation kernel is a dual of spectral_dtb function, gives 
% 'fake' spectral solutions to bounded drift diffusion and can handle arbitrary 
% changing bounds.
%
% where (all inputs in SI units)
%   'drift' is the vector of drift rates,
%   't' is the time series in seconds,
%   'Bup & Blo' are ROW vector bounds with Blo(t) < Bup(t), can be scalars if 
%       bounds are flat (+/-Inf bounds are allowed),
%   'y0' is initial distribution position,
%   'sigma' is the standard deviation,  
%   'rngseed' is the random number generator seed,
%   'notabs_flag' is flag to calculate probability of not absorped drift
%       above zero,
%   'useGPU' is a flag to use GPU for calculation or not, by default, false,
%   'trials' is the number of trials to compute.
%
%   and
%
%   D is a structure including all the above input arguments and
%   'D.up.p(drift)' is the total probability of hitting the upper bound for 
%       each drift level,
%   'D.up.mean_t(drift)' is the mean decision time for upper bound,
%   'D.up.pdf_t(t,drift)' is the probability of upper bound hit at each time 
%       (sums to D.up.p),
%   'D.up.cdf_t(t,drift)' is the cumulative probability of upper bound hit 
%       at each time (ends at D.up.p), 
%   and 
%   the '.lo' vesions as well,
%   'D.drifts' is the drifts used,
%   'D.bounds' is the bounds used,
%   'D.t' is the time used,
%   'D.useGPU' is the flag to indicate using GPU for computation or not,
%   'D.trials' is the trials for each drift.
%

% Copyright 2014 Jian Wang

if nargin < 7 || ~exist('useGPU','var')
    useGPU = false;
end

if nargin < 8 || ~exist('trials','var')
    trials = 2000;
end

% In theory, Monte Carlo simulation is developed basing on the assumption
% that the standard deviation of drift rate is 1 in 1 millisecond. Thus all
% the fitting parameters are corresponding value when time is in unit of
% second, and need to be converted for time in unit of millisecond.
drift = drift / sqrt(1000);
Bup = Bup * sqrt(1000);
Blo = Blo * sqrt(1000);
y0 = y0 * sqrt(1000);
sigma = sigma / sqrt(1000);

nt = length(t);
dt = (t(2) - t(1))*1000; % millisecond.
nd = length(drift);

% set random number generator seed
if useGPU
    parallel.gpu.rng(rngseed);
else
    rng(rngseed); 
    % For CPU, pre-generate random number matrix for reuse
    persistent rndcloud; %#ok<TLEV>
    if isempty(rndcloud)
        rndcloud = randn(nt-1,trials,nd);
    end                
end

% Expand flat bounds
if numel(Bup)==1 
    Bup = repmat(Bup,nt,1);
end

if numel(Blo)==1 
    Blo = repmat(Blo,nt,1);
end

D = struct('drift',drift,...
    't',t,...
    'Bup',Bup,...
    'Blo',Blo,...
    'y0',y0,...
    'sigma',sigma,...
    'rngseed',rngseed,...
    'useGPU',useGPU,...
    'trials',trials);

% The standard deviation in 1 millisecond is assumed to 1 technically. If
% the time step size is less than 1 millisecond, the standard deviation for
% each time step has to be multipled by sqrt(dt) so that the total variance
% for all the drifts in 1 millisecond will be equal to the standard
% deviation of 1 single one-millisecond long drift.
drift = drift * dt;
sigma = sigma * sqrt(dt);

% Monte Carlo Simulation Kernel
if useGPU
    accelerator = 'GPU';
else
    accelerator = 'CPU';
end

switch accelerator
    case 'GPU'
        k1 = parallel.gpu.CUDAKernel('find1HitBnd.ptx','find1HitBnd.cu');
        k1.ThreadBlockSize = [512, 1, 1];
        blkSizeX = int32(ceil(trials / k1.ThreadBlockSize(1)));        
        k1.GridSize = [blkSizeX, 1, 1];
        
        % _gpu subfix implies array allocated on GPU memory
        up_pdf_t_gpu = gpuArray.zeros(nt,nd);
        lo_pdf_t_gpu = gpuArray.zeros(nt,nd);
                
        BupMtx_gpu = repmat(gpuArray(Bup'),1,trials);                    
        BloMtx_gpu = repmat(gpuArray(Blo'),1,trials);
        
        for n1 = 1:nd % drift rate            
            dftForce_gpu = gpuArray.randn(nt-1,trials) * sigma(n1) + drift(n1);                                               
            dftSum_gpu = [gpuArray.zeros(1,trials); cumsum(dftForce_gpu,1)] + y0;
                                                
            dftBup_gpu = dftSum_gpu >= BupMtx_gpu;                        
            dftBlo_gpu = dftSum_gpu <= BloMtx_gpu;
                        
            clear dftSum_gpu dftForce_gpu
            hitUp_gpu = gpuArray.zeros(nt,trials);                        
            hitLo_gpu = gpuArray.zeros(nt,trials);
            
            [hitUp_gpu, hitLo_gpu] = feval(k1,dftBup_gpu,hitUp_gpu,...
                dftBlo_gpu,hitLo_gpu,int32(nt),int32(trials));
            
            hitUp_gpu = cumsum(hitUp_gpu,2);
            hitLo_gpu = cumsum(hitLo_gpu,2);
            
            up_pdf_t_gpu(:,n1) = hitUp_gpu(:,end) / trials;            
            lo_pdf_t_gpu(:,n1) = hitLo_gpu(:,end) / trials;                        
        end
                
        D.up.pdf_t = transpose(gather(up_pdf_t_gpu));        
        D.lo.pdf_t = transpose(gather(lo_pdf_t_gpu));
    
    case 'CPU'
        up_pdf_t = zeros(nt,nd); % pre-allocate
        lo_pdf_t = zeros(nt,nd);
        
        if notabs_flag
           pos_t = zeros(nt,nd); 
        end
        
        BupMtx = repmat(Bup',1,trials);
        BloMtx = repmat(Blo',1,trials);
                
        for n1 = 1:nd                        
            dftForce = rndcloud(:,:,n1) * sigma(n1) + drift(n1);
            dftSum = [zeros(1,trials); cumsum(dftForce,1)] + y0;
                                                                       
            dftBup = dftSum >= BupMtx;
            dftBlo = dftSum <= BloMtx;            
            
            hitUp = zeros(nt,1);
            hitLo = zeros(nt,1);
            
            for n2 = 1:trials
                iu = find(dftBup(:,n2),1,'first'); % index of hitting up bound
                il = find(dftBlo(:,n2),1,'first'); % index of hitting lo bound
                
                % Update hitting bound record.
                if ~isempty(iu) && isempty(il)
                    hitUp(iu) = hitUp(iu) + 1;
                    
                    if notabs_flag
                       pos_t(1:iu-1,n1) = pos_t(1:iu-1,n1) + (dftSum(1:iu-1,n1)>0);                    
                    end                    
                elseif isempty(iu) && ~isempty(il)
                    hitLo(il) = hitLo(il) + 1;
                    
                    if notabs_flag
                       pos_t(1:il-1,n1) = pos_t(1:il-1,n1) + (dftSum(1:il-1,n1)>0);                     
                    end
                elseif ~isempty(iu) && ~isempty(il)
                    if iu <= il
                        hitUp(iu) = hitUp(iu) + 1;
                        
                        if notabs_flag
                            pos_t(1:iu-1,n1) = pos_t(1:iu-1,n1) + (dftSum(1:iu-1,n1)>0);
                        end
                    else
                        hitLo(il) = hitLo(il) + 1;
                        
                        if notabs_flag
                            pos_t(1:il-1,n1) = pos_t(1:il-1,n1) + (dftSum(1:il-1,n1)>0);
                        end
                    end
                end                
            end
            
            up_pdf_t(:,n1) = hitUp / trials;
            lo_pdf_t(:,n1) = hitLo / trials;       
            
            if notabs_flag
               pos_t = pos_t / trials; 
            end
        end                                     
        
        D.up.pdf_t = up_pdf_t';
        D.lo.pdf_t = lo_pdf_t';
        
        if notabs_flag
           D.notabs.pos_t = pos_t'; 
        end
end

D.up.p = sum(D.up.pdf_t,2);
D.lo.p = sum(D.lo.pdf_t,2);    

D.up.mean_t = transpose(t'*D.up.pdf_t')./D.up.p;
D.lo.mean_t = transpose(t'*D.lo.pdf_t')./D.lo.p;

D.up.cdf_t = cumsum(D.up.pdf_t,2);
D.lo.cdf_t = cumsum(D.lo.pdf_t,2);


