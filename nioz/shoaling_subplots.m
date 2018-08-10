clear all
close all

%% Ajustable params
whichFronts = 1:4 % from 1 to 8  <------------------ Ajust here no. subplots
                  
%whichFronts = 5:9 % from 1 to 8  <------------------ Ajust here no. subplots
noFronts = length(whichFronts); 

timeWindow = 6; % hours
letterLabel = ['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'g'; 'h'; 'i'];
FS = 14;
FS1 = 12;
FS2 = 16;
figHeight = 25; %cm
%figHeight = 20; %cm
figWidth = 16; %cm
totalDepth = 919;
V = 5:.05:10;


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

%% Compute GammaVec
dt = 40/60/24;
timeVec = round(timeSensor1(1)*24)/24:dt/4:timeSensor1(end);
GammaVec = nan(size(timeVec));

for i = 1:length(timeVec)   
    I = find(timeSensor1>=timeVec(i)-dt/2 & timeSensor1<timeVec(i)+dt/2);
    GammaVec(i) = nanmean(JbVec(I))./nanmean(epsVec(I));
end
time2_gamma = time2maxV('./roc12.mat',timeVec); 


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
rigs = 0.08; % very right of figure
tops = 0.0; % top of figure
bots = 0.07; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

% $$$ % If using espwrite
% $$$ % *********************** Adjust_space.m ************************ %
% $$$ % Fields required by the function adjust_space.m. Please fill every
% $$$ % of the following and call "adjust_space" in the script whenever
% $$$ % you want. Do not touch four last fields
% $$$ ncol = 1; % no. subplot column
% $$$ nrow = noFronts; % no. subplot row
% $$$ dx = 0.01 ; % horiz. space between subplots
% $$$ dy = 0.06; % vert. space between subplots
% $$$ lefs = 0.1; % very left of figure
% $$$ rigs = 0.1; % very right of figure
% $$$ tops = 0.02; % top of figure
% $$$ bots = 0.02; % bottom of figure
% $$$ figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
% $$$ figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
% $$$ count_col = 1;
% $$$ count_row = 1;
% $$$ % *************************************************************** %

figure(1)
clf
set(gcf, 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 figWidth figHeight])

for j = 1:noFronts
    
    i = whichFronts(j);
    
    %    I = (frontsIndex(i)-2*fracWindow):(frontsIndex(i)+3*fracWindow); % <-------------------- This also could be adjusted
    I = (frontsIndex(i)-2.5*fracWindow):(frontsIndex(i)+2.5*fracWindow); % <-------------------- This also could be adjusted
    timeLims = [timeSensor1(I(2)) timeSensor1(I(end))];
    
    s = subplot(noFronts,1,j);
    imagesc(time2(I), totalDepth-zVecReg, Titp(:,I))
    hold on
    %contour(timeSensor1, totalDepth-zVecReg, Titp, V, 'color',
    %'k')
    plot([0 0], [min(totalDepth-zVecReg) max(totalDepth-zVecReg)], ...
         'color', [1 1 1], 'lineStyle', '--', 'linewidth', 2)
    hold off
    %c = colorbar
    %ti = ylabel(c, '\theta_0(^{\circ}C)', 'FontSize', FS, 'fontWeight','bold'); 
    %ylabel('hab (m)')
    %caxis([1032 1032.1]-1000)
    
    set(gca, 'xtick', -.3:.05:.3)
    if j == noFronts & whichFronts(1) ~= 1
        xlabel('\phi (days)', 'fontSize', FS, 'fontWeight', 'bold')        
    elseif j == 4 & whichFronts(1) == 1
        
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
    
    % letter ID
    text(-.12, 110, letterLabel(i), 'fontSize', 18, 'fontWeight', ...
         'bold', 'color', [1 1 1])
    
    %% Add second axis system
    II = find(timeVec>=timeLims(1) & timeVec<=timeLims(2));
    hAxes = gca;
    hAxes_pos = get(hAxes,'Position');
    hAxes_pos(4) = hAxes_pos(4)/1;
    hAxes2 = axes('Position',hAxes_pos);
    plot(time2_gamma(II), GammaVec(II),'LineWidth',2,'Color',[1 1 1]*.8)
    set(hAxes2,'YAxisLocation','right', 'Color','none', 'XTickLabel',[],'YTickLabel',[])
    h1_xlim = get(hAxes,'XLim'); % store x-axis limits of first axes
    set(hAxes2,'XLim',h1_xlim) % specify x-axis limits of second axes
    set(hAxes2,'YLim',[0 1]) % specify y-axis limits of second axes
    set(hAxes2,'Ytick',[0 .2 .4 .6 .8 1]) % specify y-axis limits of second axes
    set(hAxes2,'ycolor',[1 1 1]*.8)
    ti_y = ylabel(hAxes2, '\gamma', 'fontSize', 16, 'fontWeight', 'bold');
    set(hAxes2, 'box', 'off')
    text(h1_xlim(2)-.02, .2, '0.2', 'color', [1 1 1]*.8, 'fontSize', FS, 'fontWeight', 'bold')
    text(h1_xlim(2)-.02, .4, '0.4', 'color', [1 1 1]*.8, 'fontSize', FS, 'fontWeight', 'bold')
    text(h1_xlim(2)-.02, .6, '0.6', 'color', [1 1 1]*.8, 'fontSize', FS, 'fontWeight', 'bold')
    text(h1_xlim(2)-.02, .8, '0.8', 'color', [1 1 1]*.8, 'fontSize', FS, 'fontWeight', 'bold')
    ti_y_pos = get(ti_y, 'pos');
    ti_y_pos(1) = ti_y_pos(1)-.02;
    ti_y_pos(2) = ti_y_pos(2)+.5;
    set(ti_y, 'pos', ti_y_pos, 'rotation', 0)
    
    
    % add Gamma*
    %text(-.06, 15, sprintf('<\\gamma> = %1.2f',nanmean(JbVec(I)./epsVec(I))), 'fontSize', FS2, 'fontWeight', 'bold')
    text(-.123, .1, sprintf('<\\gamma> = %1.2f',nanmean(JbVec(I)./epsVec(I))), 'fontSize', FS2, 'fontWeight', 'bold','Color',[1 1 1]*.0)
    
end


% $$$ outfile = 'shoaling_toberename_epsutils.eps';
% $$$ epswrite(outfile)

% $$$ %outfile = sprintf('Reorder2D_Jb_20121009_%02dh25.png',myTimes(ii));
% $$$ outfile = 'shoaling_toberename.png';
% $$$ set(gcf, 'renderer', 'painters')
% $$$ print('-dpng', '-r300', outfile)
% $$$ 
% $$$ 
%outfile = sprintf('Jb-eps_%02d.eps',ii)
set(gcf, 'renderer', 'painters')
outfile = 'shoaling_toberename.eps';
print('-depsc', outfile)

