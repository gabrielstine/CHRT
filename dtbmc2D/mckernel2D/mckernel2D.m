function D =  mckernel2D(drift,t,Bup,Blo,y0,sigma,rngseed,notabs_flag,useGPU,Roh,trials)
%MCKERNEL2D Monte Carlo simulation kernel for MCFIT fitting
%
% D =  mckernel2D(drift,t,Bup,Blo,y0,sigma,rngseed,notabs_flag,useGPU,Roh,trials)
%
% This kernel applies Monte Carlo method to simulate the accumulation
% of decision information in each trial and calculate the probability
% distribution as spectral_dtbAA & FP4 methods. 
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
%   'Roh' is the (1,1) covariance in a 2 x 2 covariance matrix,
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
%
%   'D.notabs.pos_t(drifts,t)' is the probability of not being absorbed at 
%       (y>0),
%
%   'D.drifts' is the drifts used,
%   'D.bounds' is the bounds used,
%   'D.t' is the time used,
%   'D.useGPU' is the flag to indicate using GPU for computation or not,
%   'D.trials' is the trials for each drift.
%

% Copyright 2015 Jian Wang

% This Monte Carlo simulation is developed basing on the assumption that
% the standard deviation of drift rate is 1 in 1 millisecond and time is in
% the unit of 1 millisecond. Since all the fitting parameters are  
% corresponding value when time is in unit of second, all parameters need
% to be converted when time is in unit of millisecond.
drift = drift / sqrt(1000);
Bup = Bup * sqrt(1000);
Blo = Blo * sqrt(1000);
y0 = y0 * sqrt(1000);
sigma = sigma / sqrt(1000);

nt = length(t);
dt = (t(2) - t(1)) * 1000; % millisecond.
nd = length(drift);

if nargin < 7 || ~exist('useGPU','var')
    useGPU = false;
end

if nargin < 9 || ~exist('trials','var')
    trials = nt;
end

% Expand flat bounds.
if numel(Bup)==1 
    Bup = repmat(Bup,nt,1);
end

if numel(Blo)==1 
    Blo = repmat(Blo,nt,1);
end

% Set random number generator seed.
if useGPU
    parallel.gpu.rng(rngseed);
else
    rng(rngseed); 
    % Generate random number matrix for reuse.
    persistent rndPersistent; %#ok<TLEV>
    if isempty(rndPersistent)
        rndPersistent = randn(2,nd*trials*(nt-1));
    end                
end

D = struct('drift',drift * sqrt(1000),...
    't',t,...
    'Bup',Bup / sqrt(1000),...
    'Blo',Blo / sqrt(1000),...
    'y0',y0 / sqrt(1000),...
    'sigma',sigma * sqrt(1000),...
    'rngseed',rngseed,...
    'useGPU',useGPU,...
    'trials',trials,...
    'Roh',Roh);

% The standard deviation in 1 millisecond is assumed to 1 technically. If
% the time step size is less than 1 millisecond, the standard deviation for
% each time step has to be multipled by sqrt(dt) so that the total variance
% for all the drifts in 1 millisecond will be equal to the standard
% deviation of 1 single one-millisecond-long drift.
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
        if ~notabs_flag  
           % Build GPU kernel to find only bound hit.
           k1 = parallel.gpu.CUDAKernel('findBndHit2D.ptx','findBndHit2D.cu','withoutPos');                  
        else     
           % Build GPU kernel to find both bound hit and postive pdf.
           k1 = parallel.gpu.CUDAKernel('findBndHit2D.ptx','findBndHit2D.cu','withPos');
        end
        
        nthreads = 2560; % Optimized for Nvidia K20c. 
        k1.ThreadBlockSize = [512, 1, 1];                
        k1.GridSize = [nthreads / k1.ThreadBlockSize(1), 1, 1];
        n2rp = ceil(trials / nthreads);
               
        % _gpu subfix implies array allocated on GPU memory.
        up_pdf_t_gpu = gpuArray.zeros(nt,nd);
        lo_pdf_t_gpu = gpuArray.zeros(nt,nd);
        
        if notabs_flag
           pos_t_gpu = gpuArray.zeros(nt,nd); 
        end
        
        Bup_gpu = gpuArray(Bup');                    
        Blo_gpu = gpuArray(Blo');
                 
        T_cov = cholcov([1.0, Roh; Roh, 1.0]);
                        
        for n1 = 1:nd                      
            for n2 = 1:n2rp
                randpool_gpu = transpose(gpuArray.randn(nt*nthreads, 2) * T_cov * sigma(n1));                
                hitUp_gpu = gpuArray.zeros(nt,nthreads);
                hitLo_gpu = gpuArray.zeros(nt,nthreads);
                
                if notabs_flag
                    pos_gpu = gpuArray.zeros(nt,nthreads);
                end
                
                if ~notabs_flag
                    [hitUp_gpu, hitLo_gpu] = feval(k1,randpool_gpu,...
                        Bup_gpu,Blo_gpu,hitUp_gpu,hitLo_gpu,...
                        int32(nt),int32(nthreads),y0,drift(n1),Roh);
                else
                    [hitUp_gpu, hitLo_gpu, pos_gpu] = feval(k1,randpool_gpu,...
                        Bup_gpu,Blo_gpu,hitUp_gpu,hitLo_gpu,...
                        int32(nt),int32(nthreads),y0,drift(n1),Roh,pos_gpu);
                end
                
                up_pdf_t_gpu(:,n1) = up_pdf_t_gpu(:,n1) + sum(hitUp_gpu,2);
                lo_pdf_t_gpu(:,n1) = lo_pdf_t_gpu(:,n1) + sum(hitLo_gpu,2);
                
                if notabs_flag
                    pos_gpu = cumsum(pos_gpu,2);
                    pos_t_gpu(:,n1) = pos_t_gpu(:,n1) + pos_gpu(:,end);
                end
                
                clear hitUp_gpu hitLo_gpu pos_gpu;
            end
        end
        
        allThreads = nthreads * n2rp;
        
        up_pdf_t_gpu = up_pdf_t_gpu / allThreads;
        lo_pdf_t_gpu = lo_pdf_t_gpu / allThreads;      
        D.up.pdf_t = transpose(gather(up_pdf_t_gpu));        
        D.lo.pdf_t = transpose(gather(lo_pdf_t_gpu));
        
        if notabs_flag
            pos_t_gpu = pos_t_gpu / allThreads;
            D.notabs.pos_t = transpose(gather(pos_t_gpu));
        end
    
    case 'CPU'
        up_pdf_t = zeros(nt,nd);
        lo_pdf_t = zeros(nt,nd);
        
        if notabs_flag
           pos_t = zeros(nt,nd); 
        end     
        
        T_cov = cholcov([1.0, Roh; Roh, 1.0]);
        %rndTemp = transpose(rndPersistent * T_cov);
                
        for n1 = 1:nd                                   
            hitUp = zeros(nt,1);
            hitLo = zeros(nt,1);
            
            for n2 = 1:trials                
                cu = 0;
                cl = 0;
                
                for n3 = 1:nt-1
                   % Cumulative drift.
                   ldx = n3 + (n2-1)*(nt-1) + (n1-1)*(nt-1)*trials;
                   cu = cu + ...
                       (rndPersistent(1,ldx)*T_cov(1,1)+rndPersistent(2,ldx)*T_cov(2,1))...
                       *sigma(n1) + drift(n1);
                   cl = cl + ...
                       (rndPersistent(1,ldx)*T_cov(1,2)+rndPersistent(2,ldx)*T_cov(2,2))...
                       *sigma(n1) - drift(n1);
                   
                   % Check hitting lower reflactive bound.
                   if cu < Blo(1+n3)
                       cu = Blo(1+n3);
                   end
                   
                   if cl < Blo(1+n3)
                       cl = Blo(1+n3);
                   end
                   
                   % Check positive non-absorptive possibility.
                   if notabs_flag
                       if cu > cl
                           pos_t(1+n3,n1) = pos_t(1+n3,n1) +1;
                       end
                   end
                   
                   % Check & update upper bound hitting record.
                   if cu > Bup(1+n3)
                       hitUp(1+n3) = hitUp(1+n3) +1;
                       break;
                   elseif cl > Bup(1+n3)
                       hitLo(1+n3) = hitLo(1+n3) +1;
                       break;
                   end                   
                end                     
            end
            
            up_pdf_t(:,n1) = hitUp / trials;
            lo_pdf_t(:,n1) = hitLo / trials;       
        end
        
        if notabs_flag
            pos_t = pos_t / trials;
        end
                                                     
        D.up.pdf_t = up_pdf_t';
        D.lo.pdf_t = lo_pdf_t';
        
        if notabs_flag
           D.notabs.pos_t = pos_t'; 
        end
end

D.up.p = sum(D.up.pdf_t,2);
D.lo.p = sum(D.lo.pdf_t,2);    

D.up.mean_t = transpose(t'*D.up.pdf_t') ./ D.up.p;
D.lo.mean_t = transpose(t'*D.lo.pdf_t') ./ D.lo.p;

D.up.cdf_t = cumsum(D.up.pdf_t,2);
D.lo.cdf_t = cumsum(D.lo.pdf_t,2);


