clear all
close all
FS = 11;
FS2 = 13;
FS3 = 16;
% generated (manually) in APEF_vs_Thorpe.m, i.e.:
% save gammaInfo.mat timeSensor1 JbVec JbVec1 JbVec2 epsVec epsVec1 epsVec2 Reb
load gammaInfo2.mat

%%  -------------- AVERAGE GAMMA ------------- 

dt = 40/60/24;
timeVec = round(timeSensor1(1)*24)/24:dt/2:timeSensor1(end);
GammaVec = nan(size(timeVec));

for i = 1:length(timeVec)   
    I = find(timeSensor1>=timeVec(i)-dt/2 & timeSensor1<timeVec(i)+dt/2);
    GammaVec(i) = nanmean(JbVec(I))./nanmean(epsVec(I));
end

figure(2)
clf
set(gcf,'PaperUnits','centimeters', 'paperposition', [0 0 16 10])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.05; % vert. space between subplots
lefs = 0.12; % very left of figure
rigs = 0.05; % very right of figure
tops = 0.05; % top of figure
bots = 0.15; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

%% get timing of fronts
timeWindow = 6; % hours
noPtsWindow = round((timeWindow/24)./dt); % no. of pts in the window
time2 = time2maxV('./roc12.mat',timeVec)*24; 
frontsIndex =[];
for i = 2:length(time2)-1
    if abs(time2(i)) < abs(time2(i-1)) & abs(time2(i)) < abs(time2(i+1))
        frontsIndex = [frontsIndex; i];
    end
end

%% plot
plot(timeVec, GammaVec, 'k', 'linewidth', 2)
hold on
for i = 1:length(frontsIndex)
    I = (frontsIndex(i)-noPtsWindow):(frontsIndex(i)+noPtsWindow); % <-------------------- This also could be adjusted
    
    X = [timeVec(I(1)) timeVec(I(end)) timeVec(I(end)) timeVec(I(1))  timeVec(I(1))];
    Y = [0.01 0.01 .99 .99 .01];
    patch(X, Y, [1 1 1]*.8, 'linestyle', 'none');
end
plot(timeVec, GammaVec, 'k', 'linewidth', 2)
ylabel('\gamma')
datetick('x', 7)
xlabel('Oct. 2012')
set(gcf, 'renderer', 'painters')
xlim([timeVec(1) timeVec(end)])
print('-depsc2', 'GammaTimeseries.eps')


%%  -------------- HISTOGRAM ------------- 
figure(1)
clf
set(gcf,'PaperUnits','centimeters', 'paperposition', [0 0 14 18])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 3; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.05; % vert. space between subplots
lefs = 0.12; % very left of figure
rigs = 0.05; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

classEdge = [0:.05:.6];
dclass = classEdge(2)-classEdge(1);


%% S1 %%
s1 = subplot(311);
plot(timeVec, GammaVec, 'k', 'linewidth', 2)
hold on
for i = 1:length(frontsIndex)
    I = (frontsIndex(i)-noPtsWindow):(frontsIndex(i)+noPtsWindow); % <-------------------- This also could be adjusted
    
    X = [timeVec(I(1)) timeVec(I(end)) timeVec(I(end)) timeVec(I(1))  timeVec(I(1))];
    Y = [0.01 0.01 .99 .99 .01];
    patch(X, Y, [1 1 1]*.8, 'linestyle', 'none');
end
plot(timeVec, GammaVec, 'k', 'linewidth', 2)
ylabel('\gamma', 'fontSize', FS, 'fontWeight', 'bold')
datetick('x', 7)
title('Oct. 2012', 'fontSize', FS, 'fontWeight', 'bold')
xlim([timeVec(1) timeVec(end)])
set(gca, 'fontSize', FS2)
set(gca, 'ytick', [0 .25 .5 .75 1])
set(gca, 'tickdir', 'out')
adjust_space


%% S2 %%
s2 = subplot(312);
I = find(time2>=-timeWindow/2 & time2<=timeWindow/2);
EVec = GammaVec(I); % sub-sample

I = find(EVec>max(classEdge));
EVec(I) = max(classEdge);
    
clear relativeNo;
for j = 1:length(classEdge)-1
    I = find(EVec>=classEdge(j) & EVec<=classEdge(j+1));
    relativeNo(j) = length(I)/length(EVec);
end
  
bar(classEdge(1:end-1)+dclass/2, relativeNo, 'facecolor', [.2 .2 .2]);
hold on
M = nanmean(EVec);
plot([M M], [0 .22], '--', 'color', [1 1 1]*0, 'linewidth', 1)
plot([.2 .2], [0 .55], '--', 'color', [1 1 1]*0, 'linewidth', 2)
hold off
ylim([0 .3]);
xlim([classEdge(1) classEdge(end)]);
set(gca, 'xtick', classEdge(2:end-1))
set(gca, 'ytick', [0:.1:.5])
y0 = .28;
text(.58, y0, '\phi \in [-0.125, 0.125] day', 'HorizontalAlignment', 'right', 'fontSize', FS, 'fontWeight', 'bold')
text(.58, y0-.03, sprintf('<\\gamma> = %3.2f', M), 'HorizontalAlignment', 'right', 'fontSize', FS, 'fontWeight', 'bold')
text(.58, y0-.075, sprintf('\\mu_{1/2} = %3.2f', median(EVec)),'HorizontalAlignment', 'right', 'fontSize', FS, 'fontWeight', 'bold')
text(.58, y0-.1, sprintf('n = %d', length(EVec)),'HorizontalAlignment',  'right', 'fontSize', FS, 'fontWeight', 'bold')
set(gca, 'xtick', [0:.1:5])
set(gca, 'tickdir', 'out')
set(gca, 'xticklabel', '')
ylabel('Relative Freq.', 'fontSize', FS2, 'fontWeight', 'bold')
set(gca, 'fontSize', FS2)
adjust_space

%% bootstrap
nboot = 1000;
N = length(EVec);
gam_boot_b = nan(1,nboot);
for iboot = 1:nboot
    % Log-Normal distrib. bootstrap
    r = rand(N,1);
    r = ceil(r*N/1);    
    gam_boot_b(iboot) = nanmean(EVec(r));
end
CI_2p5 = round(2.5/100*nboot);
CI_97p5 = round(97.5/100*nboot);
gamSort = sort(gam_boot_b);
gam2p5 = gamSort(CI_2p5);
gam97p5 = gamSort(CI_97p5); 
[nanmean(gamSort) gam2p5 gam97p5]


%% S3 %%
s3 = subplot(313);
I = find(time2<-timeWindow/2 | time2>timeWindow/2);
EVec = GammaVec(I); % sub-sample

I = find(EVec>max(classEdge));
EVec(I) = max(classEdge);
    
clear relativeNo;
for j = 1:length(classEdge)-1
    I = find(EVec>=classEdge(j) & EVec<=classEdge(j+1));
    relativeNo(j) = length(I)/length(EVec);
end
  
bar(classEdge(1:end-1)+dclass/2, relativeNo, 'facecolor', [.2 .2 .2]);
hold on
M = nanmean(EVec);
plot([M M], [0 .22], '--', 'color', [1 1 1]*0, 'linewidth', 1)
plot([.2 .2], [0 .55], '--', 'color', [1 1 1]*0, 'linewidth', 2)
hold off
ylim([0 .32]);
xlim([classEdge(1) classEdge(end)]);
set(gca, 'xtick', classEdge(2:end-1))
set(gca, 'ytick', [0:.1:.5])

text(.58, y0, '|\phi| > 0.125 day', 'HorizontalAlignment', 'right', 'fontSize', FS, 'fontWeight', 'bold')
text(.58, y0-.03, sprintf('<\\gamma> = %3.2f', M), 'HorizontalAlignment', 'right', 'fontSize', FS, 'fontWeight', 'bold')
text(.58, y0-.07, sprintf('\\mu_{1/2} = %3.2f', median(EVec)),'HorizontalAlignment', 'right', 'fontSize', FS, 'fontWeight', 'bold')
text(.58, y0-.095, sprintf('n = %d', length(EVec)),'HorizontalAlignment',  'right', 'fontSize', FS, 'fontWeight', 'bold')

set(gca, 'xtick', [0:.1:5])
set(gca, 'tickdir', 'out')
%set(gca, 'xticklabel', '')
ylabel('Relative Freq.', 'fontSize', FS2, 'fontWeight', 'bold')
xlabel('\gamma', 'fontSize', FS2, 'fontWeight', 'bold')
set(gca, 'fontSize', FS2)
adjust_space

%% bootstrap
nboot = 1000;
N = length(EVec);
gam_boot_b = nan(1,nboot);
for iboot = 1:nboot
    % Log-Normal distrib. bootstrap
    r = rand(N,1);
    r = ceil(r*N/1);    
    gam_boot_b(iboot) = nanmean(EVec(r));
end
CI_2p5 = round(2.5/100*nboot);
CI_97p5 = round(97.5/100*nboot);
gamSort = sort(gam_boot_b);
gam2p5 = gamSort(CI_2p5);
gam97p5 = gamSort(CI_97p5); 
[nanmean(gamSort) gam2p5 gam97p5]


print('-depsc2', 'GammaHistogram2.eps')

