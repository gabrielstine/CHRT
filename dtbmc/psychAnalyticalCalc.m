function [t1,p] = psychAnalyticalCalcTest(cohs, theta)

% theta=[k,A,t1_res,t2_res]

k      = theta(1);
A      = theta(2);
t1_res = theta(3);

if length(theta) ==4;
t2_res = theta(4);
end

for i = 1:length(cohs)
C = cohs(i);

if C < 0
    t1(i,1)=A./(k.*C).*tanh(k.*C.*A)+t2_res;
elseif C > 0
    t1(i,1)=A./(k.*C).*tanh(k.*C.*A)+t1_res;
else
    t1(i,1)=A.^2+t1_res;
    %t2(i,1)=A.^2+t2_res;
end

p(i,1)=1./(1+exp(-2.*k.*C.*A));
end