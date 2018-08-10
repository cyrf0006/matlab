clear all
close all


%% Load data from APEF_vs_Thorpe
load apef_thorpe_wholeTM.mat
S2Vec = nanmean(S2bg, 1);
N2Vec = nanmean(N2bg, 1);
GamVec = JbVec./epsVec;
RiVec = N2Vec./S2Vec;
timeVec = timeSensor;

%% Few parameters
timeWindowMinutes = 10;
windowSize = timeWindowMinutes/1440; % in days
g = 9.81;

timeAve = timeVec(1):windowSize:timeVec(end);
RiAve = nan(length(timeAve), 1);
GamAve = nan(length(timeAve), 1);

for i = 1:length(timeAve)
   I = find(timeVec>=timeAve(i)-windowSize/2 & timeVec<=timeAve(i)+windowSize/2);
   RiAve(i) = nanmean(RiVec(I));
   GamAve(i) = nanmean(GamVec(I));
end


figure(1)
clf

plot(RiAve, GamAve, '.k')
xlabel('Ri')
ylabel('\Gamma')
set(gcf, 'renderer', 'painters')
print('-depsc2', 'Gamma_Ri_scatter.eps')