function P = FP4Wrapper(adv,dfu,t,Bup,Blow,notAbs)



% Convert input unit accordingly from second to millisecond. Function
% newFit801_bnd uses standard unit: second. While function FP4,
% FP4_BoundCrossZero uses millisecond. To convert, the advective rate
% (drift rate) and diffusion rate (standard deviation) need to be divided
% by sqrt(1000) and boundary height needs to be multiplied by sqrt(1000).
% This conversion formula can be determined by method of unit in a way that
% A/KC has unit while KC*A is unit free (More info on Flat Bound Fit).

adv = adv/sqrt(1000);
dfu = dfu/sqrt(1000);
t = t*1000;
Bup = Bup*sqrt(1000);
Blow = Blow*sqrt(1000);

% Define time and mesh grid resolution.
dt = t(2) - t(1);
dx = min([0.1, abs(Bup(1))/100, abs(Blow(1))/100]);

% Define xmesh and the ghost cells for bound crossing probability.
b_margin = repmat(4*dfu', [1 2]);
B = max([abs(Bup),abs(Blow)]);
xmesh = (-B-b_margin(1)+dx : dx : B+b_margin(2)-dx)';

% Adjust dx so that the mesh has zero and is therefore symmetric.
while ~any(abs(xmesh) < eps)
    [~, I] = min(abs(xmesh));
    delta_mesh = xmesh(I);
    b_margin = b_margin + delta_mesh;
    xmesh = (-B-b_margin(1)+dx : dx : B+b_margin(2)-dx)';
end

if mod(length(xmesh),2)~=1
    error('The length of xmesh must be odd.');
end

% Define a delta function on xmesh.
delta = zeros(size(xmesh));
delta(abs(xmesh)==min(abs(xmesh))) = 1;

% Find the probability density of decision variable across time for 
% stimulated and nonstimulated trials.
% Ptba = zeros(length(t)-1, 2, length(adv));             % time  * bound * signed_coh
% Pxta = zeros(length(xmesh), length(t)-1, length(adv)); % xmesh * time  * signed_coh
% Pg0a = zeros(length(t)-1, length(adv));
P = struct;

% Run FP4 to calculate Ptba and Pxta for each coherence and stimulation 
% condition.
for i = 1 : length(adv)
    mu = adv(i);
    sigma = dfu(i);
    uinit = delta;    
    b_change = [Blow - Blow(1); Bup - Bup(1)]';    
%     [~, ~, Ptb, ~, Pxt] = FP4(xmesh, uinit, mu, sigma, b_change, b_margin, dt);
    [~, ~, Ptb, Pg0, ~] = FP4(xmesh, uinit, mu, sigma, b_change, b_margin(i,:), dt);
    
%     Ptba(:,:,i) = [local_sum(Ptb(2:end,1),round(1/dt)), local_sum(Ptb(2:end,2),round(1/dt))];
%     Pxta(:,:,i) = Pxt(:,1/dt:1/dt:end);
%     Pg0a(:,i) = Pg0;

%     P.low.pdf_t(i,:) = local_sum(Ptb(2:end,1),round(1/dt));
    P.lo.pdf_t(i,:) = Ptb(2:end,1);
%     P.up.pdf_t(i,:) = local_sum(Ptb(2:end,2),round(1/dt));
    P.up.pdf_t(i,:) = Ptb(2:end,2);
    P.notabs.pos_t(i,:) = Pg0;
            
end










