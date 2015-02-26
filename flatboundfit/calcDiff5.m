function [t1,t2,p] = calcDiff5(cohs,theta,thetaFlag)
%CALCDIFF5 calculate Reaction Time & Probability of Choice for Fokker-Plank 
%   equation with the same bound hight (A = B).
%   [t1,t2,p] = calcDiff5(cohs,theta,type)
%
%   where
%       cohs is signed coherence,
%       theta is fitting parameter vector [k,A,t1_res,t2_res,du,dcoh,dpr],
%       thetaFlag is a bool vector that shows which 'theta' parameters are 
%       provided for fitting,
%       and
%       t1 and t2 is the estimated reaction time for rightward (or correct
%       choice) and leftward choice,
%       p is the estimated probability.
%       
%   Since unit(A^2) = unit(t1_res) and (KC * A) is unit-free, when changing time
%   unit from second to millisecond, for the thetaGuess, A should be devided by 
%   sqrt(1E3), and K should be multiplied by sqrt(1E3).
%
%   See also FITDIFF5.

%   Copyright Shadlen Lab 2015.

k      = theta(1);
A      = theta(2);
t1_res = theta(3);

t2_res = NaN;
if thetaFlag(4)
    t2_res = theta(4);        
end

du = 0.0; % Drift Force bias
if thetaFlag(5)
    du = theta(5);
end

dcoh = 0.0; % Motion Strength bias
if thetaFlag(6)
    dcoh = theta(6);
end

dpr = 0.0; % Vertical Probability bias   
if thetaFlag(7)
    dpr = theta(7);
end    

t1 = NaN(length(cohs),1);
t2 = NaN(length(cohs),1);
p  = NaN(length(cohs),1);

for i = 1:length(cohs)
    C = cohs(i);
    KC = k * (C+dcoh) + du;
    
    if KC ~= 0
        t1(i,1) = A / KC * tanh(KC * A) + t1_res;
        t2(i,1) = A / KC * tanh(KC * A) + t2_res;
    else
        t1(i,1) = A^2 + t1_res;
        t2(i,1) = A^2 + t2_res;
    end
    
    p(i,1) = 1 / (1 + exp(-2 * KC * A)) + dpr;
end