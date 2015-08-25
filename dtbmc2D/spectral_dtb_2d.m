function D =  spectral_dtb_2d(drift,C,t,B,yr,y0,notabs_flag)
% D =  spectral_dtb_2d(drift,C,t,B,yr,y0,notabs_flag)
% Spectral solutions to bounded drift diffusion of two races with correlated noise i.e.  2D drift to bound.
% The winner of the race is the one that reaches the upper bound first
% There is no lower bound (although we use an elastic lower bound for simulation)
% This method can handle arbitrary changing bounds.
%
% Inputs (all in SI units)
% ~~~~~~~~~~~~~~~
% drift:     matrix of drift rates with two columns for two races - so normally [-drift drift]
% C          2x2 covariane matrix of noise
% t:         time series to simulate in seconds
% B          matrix (2 col) or vector (for flat bounds) of lower and upper bounds with upper absorbing bounds and lower elastic
% yr:        vector of  DV grid values to propagate: length must be a power of 2
%            sensible range to choose is (with dt=sampling interval)
%            [min(B(1))-max(drift)*dt-4*sqrt(dt)  max(B(2))+max(drift)*dt+4*sqrt(dt)]
% y0:        vector of initial pdf (need not sum to 1)
% notabs_flag: flag as to whether (1) or not (0) to calculate the not-absorbed pdf 
%             (can requiree a lot of memory) - default is 0 if  not specified
%
% Outputs
% ~~~~~~~~~~~~~~~
% Returns D, a structure
% D.up.p(drift,race)       total probability of hitting the upper bound for each drift and race
% D.up.mean_t(drift,race)  mean decision time
% D.up.pdf_t(drift,t,race) pdf of bound hit for each drift at each time and tace (sums across time to D.up.p)
% D.up.cdf_t(drift,t,race) cumulative probability of above
%
% D.drifts            returns the drifts used
% D.bounds            returns the bounds used
% D.t                 returns the times used
%
% D.notabs.pdf(drifts,t,y1,y2) probability of not being absorbed and being at y1 in race 1 and y2 in race 2 at time t
%
% D.stop_k(drift)  number of time point simulations computed until the non-absorbed density is below threshold
% 
% This is a beta version by Daniel Wolpert.

if nargin<7
    notabs_flag=0; % flag to detetrmine whether to store the not-absorbed pdf
end

nt=length(t); %number of time points to simulate
dt=t(2)-t(1) %sampling interval

%scale drifts and covariance appropriately by dt for discrete simulation
drift_dt=drift*dt;
Cdt=C*dt;

nd=2; %number of dimensions
ndrift=size(drift,1); %number of drifts

ny=length(yr); %size of grid
[y{1},y{2}]=ndgrid(yr); % create grid of DVs

if round(log2(ny))~=log2(ny)
    error('Length of y must be a power of 2');
end

D.bounds=B;
D.drifts=drift;
D.y=yr;

%expand any flat bounds
if size(B,1)==1 
    B=ones(nt,1)*B;
end

% create fft of 2D  Gasssian with zero mean and cov Cdt - derived from page
% 302 proposition 11.26 of https://www.math.ucdavis.edu/~hunter/book/ch11.pdf

iC=inv(Cdt);
kk=[0:ny/2 -ny/2+1:-1]';
omega=2*pi*kk/range(yr);
[X{1},X{2}] = ndgrid(omega);

F = mvnpdf([X{1}(:) X{2}(:)],0,iC);

kk=((2*pi)^(nd/2))*sqrt(det(iC));
F=F*kk;

E1 = reshape(F,ny*ones(1,nd)); % this is the final fft for zero mean gaussian with our cov matrix
%%
D.up.pdf_t=zeros(ndrift,nt,nd);
if notabs_flag
    D.notabs.pdf=zeros(ndrift,nt,ny,ny);
end

p_threshold=10e-4; %threshold for proportion un-terminated to stop simulation

for q=1:ndrift  % iterate over drifts
    q
    u=y0;     %u is the pdf an we set u to initial state
    
    %shift mean of gaussian by drift
    dd=drift_dt(q,:)';
    
    E2=E1;
    for d=1:nd
        E2=E2.*exp(-1i*X{d}*dd(d));
    end
    
    pup=zeros(nt,nd);

    D.stop_k(q)=nt; %the maximum stop step by default
    
    for k=1:nt %iterate over time
        s=sum(u(:)); %used later for normalization
        
        %fft our current pdf
        ufft=fft2(u);
        
        %convolve with gaussian with drift in the frequency domain
        %by performing pointwide multiplication
        ufft=E2.*ufft;
        
        %turn back into time domain
        u=real(ifft2(ufft));
        u=max(u,0);
        
        s1=sum(u(:));
        u=u*s/s1; %to preserve density
        
        %select density that has crossed bounds
        for d=1:nd
            s=u.*(y{d}>B(k,2));    %this is  density that has croseed upper bound   
            u=u.*(y{d}<=B(k,2));   %remove density that has croseed upper bound 
            
            s1=u.*(y{d}<B(k,1));       % density less than lower bound
            u =u.*(y{d}>=B(k,1));    %keep only density greater than lower bounds
            
            %now add negative density *that we have removed) to the lower bound so reflective
            %(or maybe elastic)
            
            ii=find(yr>=B(k,1),1); % index of first grid above lower bound

            %add density back in along the correct dimension
            if d==1
                u(ii,:)=u(ii,:)+sum(s1,1);
            elseif d==2
                u(:,ii)=u(:,ii)+sum(s1,2);
            end
            
            %terminated density
            pup(k,d)=sum(s(:));
            D.up.pdf_t(q,k,d)=pup(k,d);
        end
        
        if notabs_flag
            D.notabs.pdf(q,k,:,:)=u;
        end
        
        %exit if our threshold is reached
        if sum(u(:))<p_threshold,
            D.stop_k(q)=k;
            break
        end
    end
    
    D.up.p(q,:)=sum(pup);
    D.up.mean_t(q,:)=(t'*pup)./D.up.p(q,:);
end

D.up.cdf_t=cumsum(D.up.pdf_t,2);
D.t=t;



