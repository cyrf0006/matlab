clear all
close all

%% Ajustable params
whichFronts = 1:4 % from 1 to 8  <------------------ Ajust here no. subplots
whichFronts = 6:8 % from 1 to 8  <------------------ Ajust here no. subplots
%whichFronts = 5:9 % from 1 to 8  <------------------ Ajust here no. subplots
noFronts = length(whichFronts); 

timeWindow = 4; % hours
letterLabel = ['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'g'; 'h'; 'i'];
FS = 12;
FS1 = 10;
FS2 = 16;
figHeight = 16; %cm
%figHeight = 20; %cm
figWidth = 16; %cm
totalDepth = 919;
V = 5:.1:10;
V1 = 5:.01:10;


%% Load whole timeserie
load apef_thorpe_wholeTM.mat
dt = (timeSensor1(2) - timeSensor1(1))*86400; % sec.
noPtsWindow = (timeWindow*60*60)./dt; % no. of pts in the window
fracWindow = round(noPtsWindow./5);

%% get timing of fronts
time2 = time2maxV('./roc12.mat',timeSensor1); 
frontsIndex =[];
for i = 2:length(time2)-1
    if abs(time2(i)) < abs(time2(i-1)) & abs(time2(i)) < abs(time2(i+1))
        frontsIndex = [frontsIndex; i];
    end
end
% discard 1st 


%% Plot
% If using regular print
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = noFronts; % no. subplot row
dx = 0.01 ; % horiz. space between subplots
dy = 0.02; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.1; % very right of figure
tops = 0.02; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

figure(1)
clf
set(gcf, 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 figWidth figHeight])

for j = 1:noFronts
    
    i = whichFronts(j);
    
    I = (frontsIndex(i)-2*fracWindow):(frontsIndex(i)+3*fracWindow); % <-------------------- This also could be adjusted
    
    s = subplot(noFronts,1,j);
    contourf(time2(I), totalDepth-zVecReg, Titp(:,I), V1, 'linestyle', 'none')
    hold on
    contour(time2(I), totalDepth-zVecReg, Titp(:,I), V, 'color','k')
% $$$     plot([0 0], [min(totalDepth-zVecReg) max(totalDepth-zVecReg)], ...
% $$$          '--k', 'linewidth', 2)
    hold off
    %c = colorbar
    %ti = ylabel(c, '\theta_0(^{\circ}C)', 'FontSize', FS, 'fontWeight','bold'); 
    %ylabel('hab (m)')
    %caxis([1032 1032.1]-1000)
    set(gca, 'xtick', -.3:.05:.3)
    if j == noFronts
        xlabel('\phi (days)', 'fontSize', FS, 'fontWeight', 'bold')        
    else
        set(gca, 'xticklabel', []);        
    end

    set(gca, 'ydir', 'normal');
    set(gca, 'tickdir', 'out')
    xlim([time2(I(1)) time2(I(end))])
    ylim([0 max(totalDepth-zVecReg)])
    ylabel('hab (m)', 'fontSize', FS, 'fontWeight', 'bold')
    c = colorbar;
    adjust_space
    drawnow
    pause(.1)

    CLIM = get(c, 'ylim');
    set(c, 'ytick', 6.5:.5:9)
    %set(c, 'ytick', round(CLIM*2)/2)
    set(c, 'clim', CLIM);
    
    cbPos = get(c, 'pos');
    cbPos(1) = cbPos(1)-.02;
    cbPos(2) = cbPos(2)+.03;
    cbPos(3) = cbPos(3)*.7;
    cbPos(4) = cbPos(4)-.04;
    set(c, 'pos', cbPos)
    ti = ylabel(c, '\theta_0(^{\circ}C)', 'FontSize', FS1, 'fontWeight','bold'); 
    ti_pos = [1 CLIM(1)-diff(CLIM)/10 1]; 
    set(ti, 'rot', 0) 
    set(ti, 'pos', ti_pos)
    drawnow
    set(gca, 'fontSize', FS)
         
    % add Gamma*
% $$$     text(.068, 15, sprintf('<\\gamma> = %1.2f',nanmean(JbVec(I)./ ...
% $$$                                                       epsVec(I))), ...
% $$$          'fontSize', FS2, 'fontWeight', 'bold', 'color',[1 1 1])
    
end

outfile = 'shoaling_EM.png';
print('-dpng', '-r300', outfile)

