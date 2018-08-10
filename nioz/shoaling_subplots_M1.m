clear all
close all

%% Ajustable params
whichFronts = 1:4 % from 1 to 8  <------------------ Ajust here no. subplots
%whichFronts = 5:9 % from 1 to 8  <------------------ Ajust here no. subplots
noFronts = length(whichFronts); 

timeWindow = 4; % hours
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
t0 = datenum(2012, 10, 7, 19, 31, 0);
%t0 = datenum(2012, 10, 7, 20, 31, 0);


%% Plot
% If using regular print
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = noFronts; % no. subplot row
dx = 0.01 ; % horiz. space between subplots
dy = 0.05; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.08; % very right of figure
tops = 0.0; % top of figure
bots = 0.07; % bottom of figure
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
    
    tRef = t0+(i-1)*24.84/24;    


    [Y, IRef] = min(abs(timeSensor1-tRef));
    I = (IRef-2*fracWindow):(IRef+3*fracWindow); % <-------------------- This also could be adjusted
    
    s = subplot(noFronts,1,j);
    imagesc(timeSensor1(I), totalDepth-zVecReg, Titp(:,I))
    hold on
    plot([tRef tRef], [min(totalDepth-zVecReg) max(totalDepth-zVecReg)], ...
         'color', [1 1 1], 'lineStyle', '--', 'linewidth', 2)
    hold off
    datetick('x', 13)
    if j == noFronts & whichFronts(1) ~= 1
        xlabel('\phi (days)', 'fontSize', FS, 'fontWeight', 'bold')        
    else
        %set(gca, 'xticklabel', []);        
    end

    set(gca, 'ydir', 'normal');
    set(gca, 'tickdir', 'out')
    xlim([timeSensor1(I(1)) timeSensor1(I(end))])
    ylim([0 max(totalDepth-zVecReg)])
    ylabel('hab (m)', 'fontSize', FS, 'fontWeight', 'bold')
    c = colorbar;
    caxis([6.8 8.2])
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
    text(0.09, 15, letterLabel(i), 'fontSize', FS2, 'fontWeight', ...
         'bold', 'color', [1 1 1])
    
    % add Gamma*
    text(-.06, 15, sprintf('<\\gamma> = %1.2f',nanmean(JbVec(I)./epsVec(I))), 'fontSize', FS2, 'fontWeight', 'bold')
    
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
outfile = 'shoaling_toberename.eps';
print('-depsc', outfile)

