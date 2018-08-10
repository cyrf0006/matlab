V = 5:.1:10;
FS = 12;
FS2 = 16;

%% Few plots
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 4; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.12; % very left of figure
rigs = 0.14; % very right of figure
tops = 0.06; % top of figure
bots = 0.07; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %
figure(1)
clf
set(gcf, 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 16 18])
s1 = subplot(4,1,[1 2]);
imagesc(timeSensor1, totalDepth-zVecReg, Titp)
hold on
contour(timeSensor1, totalDepth-zVecReg, Titp, V, 'color', 'k')
plot([datenum(2012,10,10,21,45,0) datenum(2012,10,10,21,45,0)], [0 120], 'color', [1 1 1], 'lineStyle', '--', 'lineWidth', 2)
hold off
c = colorbar('location', 'northoutside');
ti = ylabel(c, '\theta_0(^{\circ}C)', 'FontSize', FS, 'fontWeight','bold'); 
ylabel('hab (m)', 'fontSize', FS,'fontWeight','bold')
datetick
set(gca, 'xticklabel', []);
set(gca, 'ydir', 'normal');
set(gca, 'tickdir', 'out')
xlim([timeSensor1(1) timeSensor1(end)])
ylim([0 max(totalDepth-zVecReg)])

%text(timeSensor1(end)-1/1440, 110, sprintf('<Re> = %1.0e',Re),'fontSize',FS2,'horizontalAlignment', 'right','verticalAlignment','bottom','color',[1 1 1]*1,'fontWeight','bold')
text(timeSensor1(end)-1/1440, 110, '<Re> = 2x10^7','fontSize',FS2,'horizontalAlignment', 'right','verticalAlignment','bottom','color',[1 1 1]*1,'fontWeight','bold')

adjust_space
pos1 = get(gca, 'pos');
adjust_space
pos2 = get(gca, 'pos');
pos2(4) = 2*pos2(4)+dy;
set(gca, 'pos', pos2)
drawnow

cbPos = get(c, 'pos');
cbPos(1) = cbPos(1)+.08;
cbPos(3) = cbPos(3)-.1;
cbPos(2) = cbPos(2)-.04;
cbPos(4) = cbPos(4)*.4;
set(c, 'pos', cbPos)
set(ti, 'rotation', 0)
CLIM = get(gca, 'clim');
ti_pos = [CLIM(1)-diff(CLIM)/12, 1 1];
set(ti, 'pos', ti_pos)
text(timeSensor1(1)+2/1440, 115, 'a','fontSize',FS2,'fontWeight', 'bold', 'BackgroundColor',[1 1 1]);
set(gca, 'fontSize', FS)
drawnow

s2 = subplot(4,1,3);
semilogy(timeSensor1, JbVec, 'k', 'linewidth', 2) 
hold on
plot([timeSensor(1) timeSensor(end)], [1 1]*nanmean(JbVec), '--k')
plot([timeSensor(1) timeSensor(end)], [1 1]*nanmean(epsVec), 'color',  [0 1 1]*.7, 'linestyle', '--')
semilogy(timeSensor1, epsVec, 'color', [0 1 1]*.7, 'linewidth', 2)
semilogy(timeSensor1, JbVec, 'k', 'linewidth', 2) 
datetick
set(gca, 'xticklabel', []);
set(gca, 'tickdir', 'out')
xlim([timeSensor1(1) timeSensor1(end)])
ylabel('(m^2 s^{-3})', 'fontSize', FS,'fontWeight','bold')
text(timeSensor1(end)+5/1440, 1e-7, '\epsilon', 'color', [0 1 1]*.7,'fontSize',FS2,'fontWeight', 'bold')
text(timeSensor1(end)+5/1440, 6e-9, 'J^*_b', 'color', 'k', 'fontSize',FS2,'fontWeight', 'bold')
set(gca, 'ytick', [1e-9 1e-8 1e-7 1e-6])
%ylim([1e-9 5e-7])
ylim([1e-9 1e-5])
text(timeSensor1(1)+2/1440, 2e-6, 'b','fontSize',FS2,'fontWeight', 'bold', 'BackgroundColor',[1 1 1]);
adjust_space
set(gca, 'fontSize', FS)
drawnow

s3 = subplot(4,1,4);
plot(timeSensor1,JbVec./epsVec, 'k', 'linewidth', 2)
hold on
%plot(timeSensor1,JbVec*0+.2, 'color', [1 1 1]*0,'linestyle', '--')
plot(timeSensor1,JbVec*0+nanmean(JbVec./epsVec), 'color', [1 1 1]*0,'linestyle', '--')
ylim([0 1])
datetick('x', 15)
xlim([timeSensor1(1) timeSensor1(end)])
xlabel(sprintf(datestr(timeSensor1(1), 1)), 'fontSize', FS,'fontWeight','bold');
ylabel('\gamma', 'fontSize', FS,'fontWeight','bold');
text(timeSensor1(end)-42/1440, .9, sprintf('<\\gamma> = %1.2f', ...
                                           nanmean(JbVec./epsVec)), 'fontSize',FS,'fontWeight','bold')
set(gca, 'tickdir', 'out')
text(timeSensor1(1)+2/1440, .85, 'c','fontSize',FS2,'fontWeight', 'bold', 'BackgroundColor',[1 1 1]);
adjust_space
set(gca, 'fontSize', FS)

outfile = 'Jb-eps_toberename.eps';
print('-depsc2', outfile)


