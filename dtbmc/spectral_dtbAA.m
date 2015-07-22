function D =  spectral_dtbAA(drift,t,Bup,Blo,y,y0,notabs_flag,useGPU,varargin)
% D =  spectral_dtbAA(drift,t,Bup,Blo,y,y0,notabs_flag,useGPU)
%
% This is the antialiased version of spectral_dtb function, gives spectral 
% solutions to bounded drift diffusion and can handle arbitrary changing bounds.
% The evidence space is antialiased by removing the 'spiking' phenomenon.
%
% Inputs: (all in SI units)
% 'drift' is the vector of drift rates,
% 't' is the time series in seconds,
% 'Bup & Blo' are ROW vector bounds with Blo(t) < Bup(t), can be scalars if bounds 
%     are flat (+/-Inf bounds are allowed),
% 'y' is the vector of values to propagate, whose length must be a power of 2.
%     The sensible range to choose is (with dt=sampling interval)
%     [min(Blo)-max(drift)*dt-4*sqrt(dt)  max(Bup)+max(drift)*dt+4*sqrt(dt)],
% 'y0' is the vector of initial pdf (need not sum to 1),
% 'notabs_flag' is the flag of whether (1) or not (0) to calculate the notabs
%     pdf. This calculation can take a lot of memory. By default 0 if not
%     specified,
% 'useGPU' is a flag to use GPU for calculation or not, by default false.
%
% Outputs:
% Returns a structure D including all the above input arguments and
% 'D.up.p(drift)' is the total probability of hitting the upper bound for each 
%     drift level,
% 'D.up.mean_t(drift)' is the mean decision time for upper bound,
% 'D.up.pdf_t(t,drift)' is the probability of upper bound hit at each time (sums
%     to D.up.p),
% 'D.up.cdf_t(t,drift)' is the cumulative probability of upper bound hit at each
%     time (ends at D.up.p),
% and the '.lo' vesions as well,
%
% 'D.drifts' is the drifts used,
% 'D.bounds' is the bounds used,
% 'D.t' is the time used,
%
% 'D.notabs.pdf(drifts,y,t)' is the probability of not being absorbed at (y,t),
% 'D.notabs.pos_t(drifts,t)' is the probability of not being absorbed at (y>0),
% 'D.notabs.neg_t(drifts,t)' is the probability of not being absorbed at (y<0),
% 'D.notabs.y' is the sets of 'y' considered for the pdf,
% 'D.useGPU' is the flag to indicate using GPU for computation or not.

% Beta version 1.0 by Daniel Wolpert.
% Vectorized by Luke Woloszyn.
% Indexing modified and antialiased by Yul Kang.
%
% Copyright 2013 Daniel Wolpert
% Copyright 2013 Luke Woloszyn
% Copyright 2013 Yul Kang
% Copyright 2013 Jian Wang

if nargin < 7 || ~exist('notabs_flag','var')
    notabs_flag = 0; % Flag to detetrmine whether to store the notabs pdf
end

if nargin < 8 || ~exist('useGPU','var')
    useGPU = false;
end

nt = length(t);
dt = t(2)-t(1); % Sampling time interval
nd = length(drift);
ny = length(y);

if length(varargin) == 2 % diffusion term
    dfu = varargin{1};
    rngseed = varargin{2};
    rng(rngseed);
    rndnum = normrnd(0,1,[nt,1,length(dfu)]); % [nt,1,length(df)]
    rndnum = repmat(rndnum,[1,ny,1]); % [nt,ny,length(df)]
    dfu = reshape(dfu,[1,1,length(dfu)]);
    dfu = repmat(dfu,[nt,ny,1]);    
    dfu = bsxfun(@times,rndnum,dfu);   
end

if round(log2(ny)) ~= log2(ny)
    error('Length of y must be a power of 2');
end

if numel(y0) ~= numel(y)
    error('Length of y must be sames as y0');
end

% Expand any flat bounds
if numel(Bup)==1 
    Bup = repmat(Bup,nt,1);
end
if numel(Blo)==1 
    Blo = repmat(Blo,nt,1);
end

if exist('dfu','var')
    D = struct('drift',drift,'t',t,'Bup',Bup,'Blo',Blo,'y',y,'y0',y0,...
        'notabs_flag',notabs_flag,'useGPU',useGPU,'diffusion',dfu);
else
    D = struct('drift',drift,'t',t,'Bup',Bup,'Blo',Blo,'y',y,'y0',y0,...
        'notabs_flag',notabs_flag,'useGPU',useGPU);
end

% Create fft of unit variance zero mean Gaussian, repmat it to batch over
% drifts (each column corresponds to different drift)
kk = repmat([0:ny/2 -ny/2+1:-1]',[1,nd]);
omega = 2*pi*kk/range(y);
% fft of the normal distribution - scaled suitably by dt
E1 = exp(-0.5*dt*omega.^2); 
% This is the set of shifted gaussians in the frequency domain with one gaussian 
% per drift (each column)
E2 = E1.*exp(-1i.*omega.*repmat(drift,[ny,1])*dt);

% FFT of the diffusion term
if exist('dfu','var')
    tempOmega = repmat(reshape(omega,[1,size(omega)]),[nt,1,1]);
    DFUfft = exp(-1i.*tempOmega.*dfu*dt); %TODO: sqrt(dt)
end

% Preallocate
D.up.pdf_t = zeros(nd,nt);
D.lo.pdf_t = zeros(nd,nt);
if notabs_flag
    D.notabs.pdf = zeros(nd,ny,nt);
end
if useGPU
    % _gpu subfix implies array allocated on GPU memory
    D_up_pdf_t_gpu = gpuArray.zeros(nd,nt); 
    D_lo_pdf_t_gpu = gpuArray.zeros(nd,nt);
    if notabs_flag
        D_notabs_pdf_gpu = gpuArray.zeros(nd,ny,nt);
    end
end

% Initial state, repmated for batching over drifts (each column will correspond
% to different drift)
U = repmat(y0, [1,nd]);

% Repmat this too
% Y = repmat(y,[1,nd]);

p_threshold = 1.0E-5; % Threshold for proportion un-terminated to stop simulation

% Prepare anti-aliasing
[ixAliasUp, wtAliasUp] = bsxClosest(Bup, y);
[ixAliasDn, wtAliasDn] = bsxClosest(Blo, y);

dy        = y(2)-y(1); % y should be uniformly spaced.
wtAliasUp = (wtAliasUp / dy) - 0.5; 
wtAliasDn = (wtAliasDn / dy) + 0.5;

if ~useGPU
   upVA = zeros(nd,nt);
   dnVA = zeros(nd,nt);
else
   wtAliasUp_gpu = gpuArray(wtAliasUp);
   wtAliasDn_gpu = gpuArray(wtAliasDn);
   E2_gpu = gpuArray(E2);
   
   if exist('DFUfft','var')
       DFUfft_gpu = gpuArray(DFUfft);
   end
   
   U_gpu = gpuArray(U);   
   upVA_gpu = gpuArray.zeros(nd,nt); % gpuArray of upV_gpu and dnV_gpu
   dnVA_gpu = gpuArray.zeros(nd,nt); 
end

% Iterate over time
for k=1:nt 
    
    %k
    
    % FFT of current pdf
    if ~useGPU
        Ufft = fft(U);
    else        
        Ufft_gpu = fft(U_gpu);
    end
    
    % Convolve with gaussian via pointwise multiplication in frequency domain
    if ~useGPU
        Ufft = E2.*Ufft;
                
        if exist('DFUfft','var') % blur with the diffusion term
           Ufft = squeeze(DFUfft(k,:,:)).* Ufft;
        end        
    else
        Ufft_gpu = E2_gpu .* Ufft_gpu;
        
        if exist('DFUfft_gpu','var')
           Ufft_gpu = squeeze(DFUfft_gpu(k,:,:)).*Ufft_gpu; 
        end
    end
    
    % Back into time domain
    if ~useGPU
        U = max(real(ifft(Ufft)),0);
    else
        U_gpu = max(real(ifft(Ufft_gpu)),0);
    end
    
    % Select density that has crossed bounds
    if ~useGPU
        D.up.pdf_t(:,k) = sum(U(ixAliasUp(k):end,:),1); 
        D.lo.pdf_t(:,k) = sum(U(1:ixAliasDn(k)  ,:),1);
    else
        D_up_pdf_t_gpu(:,k) = sum(U_gpu(ixAliasUp(k):end,:),1);        
        D_lo_pdf_t_gpu(:,k) = sum(U_gpu(1:ixAliasDn(k)  ,:),1);               
    end
            
    % On the boundary, anti-alias
    if ~useGPU
        upV = U(ixAliasUp(k),:);
        upVA(:,k) = upV;
        dnV = U(ixAliasDn(k),:);       
        dnVA(:,k) = dnV;
    else
        upV_gpu = U_gpu(ixAliasUp(k),:);
        upVA_gpu(:,k) = upV_gpu;
        dnV_gpu = U_gpu(ixAliasDn(k),:);
        dnVA_gpu(:,k) = dnV_gpu;
    end
    
    % Keep only density within bounds
    if ~useGPU
        U(y<=Blo(k) | y>=Bup(k),:) = 0;    
    else                
        thisIdx = y<=Blo(k) | y>=Bup(k);
        U_gpu(thisIdx,:) = 0;
    end    
        
    if ~useGPU
        U(ixAliasUp(k),:) = upV * (-wtAliasUp(k));
        U(ixAliasDn(k),:) = dnV *  wtAliasDn(k);    
    else
        U_gpu(ixAliasUp(k),:) = upV_gpu * (-wtAliasUp_gpu(k)); 
        U_gpu(ixAliasDn(k),:) = dnV_gpu * wtAliasDn_gpu(k);
    end
    
    % Save if requested
    if notabs_flag
        if ~useGPU
            D.notabs.pdf(:,:,k) = U';
        else
            D_notabs_pdf_gpu(:,:,k) = transpose(U_gpu);
        end
    end
    
    % Exit if our threshold is reached (because of the batching, we have to wait 
    % until all densities have been absorbed)     
    if ~useGPU
        if sum(sum(U,1)<p_threshold)==nd 
            break;
        end
    else
        if sum(sum(U_gpu,1)<p_threshold)==nd
            break;
        end
    end    
end

% Anti-aliasing continued
if ~useGPU    
    wtAliasUp = repmat(wtAliasUp,nd,1);
    D.up.pdf_t = D.up.pdf_t + wtAliasUp .* upVA;
    wtAliasDn = repmat(wtAliasDn,nd,1);
    D.lo.pdf_t = D.lo.pdf_t - wtAliasDn .* dnVA;    
else
    wtAliasUp_gpu = repmat(wtAliasUp_gpu,nd,1);
    D_up_pdf_t_gpu = D_up_pdf_t_gpu + wtAliasUp_gpu .* upVA_gpu;
    wtAliasDn_gpu = repmat(wtAliasDn_gpu,nd,1);
    D_lo_pdf_t_gpu = D_lo_pdf_t_gpu - wtAliasDn_gpu .* dnVA_gpu;    
    D.up.pdf_t = gather(D_up_pdf_t_gpu);
    D.lo.pdf_t = gather(D_lo_pdf_t_gpu);
    
    if notabs_flag
        D.notabs.pdf = gather(D_notabs_pdf_gpu);
    end    
end

if notabs_flag
    D.notabs.pos_t = squeeze(sum(D.notabs.pdf(:,y'>=0,:),2));
    D.notabs.neg_t = squeeze(sum(D.notabs.pdf(:,y'< 0,:),2)); 
end

D.up.p = sum(D.up.pdf_t,2);
D.lo.p = sum(D.lo.pdf_t,2);
    
t = t(:);

D.up.mean_t = transpose(t'*D.up.pdf_t')./D.up.p;
D.lo.mean_t = transpose(t'*D.lo.pdf_t')./D.lo.p;

D.up.cdf_t = cumsum(D.up.pdf_t,2);
D.lo.cdf_t = cumsum(D.lo.pdf_t,2);


