function P = FP4Wrapper(adv,dfu,t,Bup,Blo,y0,notAbs_flag)
%P = FP4Wrapper(adv,dfu,t,Bup,Blow,y0,notAbs_flag)
%
%   FP4 function wrapper.
%

P = struct('drift',adv,'t',t,'Bup',Bup,'Blo',Blo);

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
Blo = Blo*sqrt(1000);
y0 = y0/sqrt(1E3);

% Define time and mesh grid resolution.
dt = t(2) - t(1);
dx = min([0.1, abs(Bup(1))/100, abs(Blo(1))/100]);

% Define xmesh and the ghost cells for bound crossing probability.
b_margin = repmat(4*dfu', [1 2]);
B = max([abs(Bup),abs(Blo)]);
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
i1 = find(xmesh>=y0, 1,'first');
i2 = find(xmesh<=y0, 1,'last');

if i1 == i2
    delta(i1) = 1;
else    
    w1 = abs(xmesh(i2)-y0);    
    w2 = abs(xmesh(i1)-y0);
    w1 = w1/(w1+w2);
    w2 = (1-w1);
    delta(i1) = w1;
    delta(i2) = w2;
end

for i = 1 : length(adv)
    mu = adv(i);
    sigma = dfu(i);
    uinit = delta;    
    b_change = [Blo - Blo(1); Bup - Bup(1)]';    
    [~, ~, Ptb, Pg0, ~] = FP4(xmesh, uinit, mu, sigma, b_change, b_margin(i,:), dt);    
    P.lo.pdf_t(i,:) = Ptb(2:end,1);
    P.up.pdf_t(i,:) = Ptb(2:end,2);
    P.notabs.pos_t(i,:) = Pg0;            
end










