% test_confMap.m

% Define the set of coherences and their relative sampling frequencies.
% This is effectively the prior over difficulty

% choose coherence set. Thesea are drift rates really (so k*c).
c = [0;.5];
c = [-1:.1:1]';
c = [0 logspace(log10(.032),log10(.512),5)]'
% c = linspace(.1,2,6)'


% weights on the coherences.
% w = normpdf(c,0,.5)
w = ones(size(c));
% w = normpdf(c,1,.5)

figure(1), clf
plot(c,w,'k-o')

% call the function
[pSet,M] = confMap(c,w)


figure(2), clf
imagesc(M.CMAP)

M.c

