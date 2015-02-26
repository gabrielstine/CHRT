function scoh = scohDynSample
% Sample code to generate random signed coherence

mu = 0.256;
sd = 0.1;

scohRange = [0.0,0.512];

scoh = -2.0;

while scoh < scohRange(1) || scoh > scohRange(2)
    scoh = normrnd(mu,sd,1,1);        
end

