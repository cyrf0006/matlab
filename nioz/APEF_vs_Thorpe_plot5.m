FS = 12;
FS2 = 16;
FS3 = 14;

%% Few plots
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 3; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.06; % vert. space between subplots
lefs = 0.23; % very left of figure
rigs = 0.16; % very right of figure
tops = 0.06; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %
figure(1)
clf
set(gcf, 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 16 14])
s1 = subplot(3,1,[1 2]);
%imagesc(timeSensor1, zVecReg, Ritp1-1000)
imagesc(timeSensor1, totalDepth-zVecReg, Titp)
hold on
contour(timeSensor1, totalDepth-zVecReg, Titp, V, 'color', 'k')
hold off
c = colorbar('location', 'northoutside');
%ylabel(c, '\sigma_1 (kg m^{-3})', 'FontSize', 10); 
ti = ylabel(c, '\theta_0(^{\circ}C)', 'FontSize', FS, 'fontWeight','bold'); 
%ylabel('hab (m)')
%caxis([1032 1032.1]-1000)
datetick
set(gca, 'xticklabel', []);
set(gca, 'yticklabel', []);
set(gca, 'ydir', 'normal');
set(gca, 'tickdir', 'out')
xlim([timeSensor1(1) timeSensor1(end)])
ylim([0 max(totalDepth-zVecReg)])

%text(timeSensor1(1)+.5/1440, 5, sprintf('<Ri> = %3.2f',nanmean(N2Vec./S2Vec)),'horizontalAlignment', 'left','verticalAlignment','bottom','color',[1 1 1]*1,'fontWeight','bold')
%text(timeSensor1(end)-1/1440, 110, sprintf('<Re> = %1.0e',Re),'fontSize',FS2,'horizontalAlignment', 'right','verticalAlignment','bottom','color',[1 1 1]*1,'fontWeight','bold')
% OVT
text(timeSensor1(end)-1/1440, 110, '<Re> = 3x10^6','fontSize',FS2,'horizontalAlignment', 'right','verticalAlignment','bottom','color',[1 1 1]*1,'fontWeight','bold')
% LOW1
% $$$ text(timeSensor1(end)-1/1440, 110, '<Re> = 8x10^6','fontSize',FS2,'horizontalAlignment', 'right','verticalAlignment','bottom','color',[1 1 1]*1,'fontWeight','bold')
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
XTICKS = get(gca, 'xtick');
drawnow

s2 = subplot(3,1,3);
semilogy(timeSensor1, JbVec, 'k', 'linewidth', 2) 
hold on
plot([timeSensor1(1) timeSensor1(end)], [1 1]*nanmean(JbVec), '--k')
plot([timeSensor1(1) timeSensor1(end)], [1 1]*nanmean(epsVec), 'color',  [0 1 1]*.7, 'linestyle', '--')
semilogy(timeSensor1, epsVec, 'color', [0 1 1]*.7, 'linewidth', 2)
semilogy(timeSensor1, JbVec, 'k', 'linewidth', 2) 

xlim([timeSensor1(1) timeSensor1(end)])
%xlabel(sprintf(datestr(timeSensor(I(1)), 1)));
ylab = ylabel('(m^2 s^{-3})', 'fontSize', FS, 'fontWeight', 'bold');
% OVT
text(timeSensor1(end)+1/1440, 1e-7, '\epsilon', 'color', [0 1 1]*.7,'fontSize',FS2,'fontWeight', 'bold')
text(timeSensor1(end)+1/1440, 6e-9, 'J^*_b', 'color', 'k','fontSize',FS3,'fontWeight', 'bold')
% LOW1
% $$$ text(timeSensor1(end)+1/1440, 1e-7, '\epsilon', 'color', [0 1 1]*.7,'fontSize',FS2,'fontWeight', 'bold')
% $$$ text(timeSensor1(end)+1/1440, 4e-9, 'J^*_b', 'color', 'k','fontSize',FS2,'fontWeight', 'bold')
set(gca, 'ytick', [1e-9 1e-8 1e-7 1e-6])
%set(gca, 'ytick', XTICKS)
%ylim([1e-9 5e-7])
ylim([1e-9 1e-5])
datetick('x', 15)
xlim([timeSensor1(1) timeSensor1(end)])
xlabel(sprintf(datestr(timeSensor1(1), 1)), 'fontSize', FS, 'fontWeight', 'bold');

set(gca, 'fontSize', FS)
ylab_pos = get(ylab, 'pos');

ylab_pos(2) = 3e-8;
ylab_pos(1) = ylab_pos(1)-.001;

set(ylab, 'pos', ylab_pos);
text(timeSensor1(1)+2/1440, 3e-6, 'b','fontSize',FS2,'fontWeight', 'bold', 'BackgroundColor',[1 1 1]);
text(timeSensor1(end)-14/1440, 3e-6, sprintf('\\gamma = %1.2f',nanmean(JbVec)./nanmean(epsVec)),'fontSize',FS2,'fontWeight','bold')
set(gca, 'tickdir', 'out')

adjust_space
drawnow

% $$$ 
% $$$ s3 = subplot(4,1,4);
% $$$ plot(timeSensor1,JbVec./epsVec, 'k', 'linewidth', 2)
% $$$ %datetick('x', 7)
% $$$ hold on
% $$$ plot(timeSensor1,JbVec*0+nanmean(JbVec./epsVec), 'color', [1 1 1]*0,'linestyle', '--')
% $$$ ylim([0 1])
% $$$ datetick('x', 15)
% $$$ xlim([timeSensor1(1) timeSensor1(end)])
% $$$ xlabel(sprintf(datestr(timeSensor1(1), 1)), 'fontSize', FS, 'fontWeight', 'bold');
% $$$ %xlabel('2012-10-11')
% $$$ %xlabel('October 2012')
% $$$ ylabel('\gamma', 'fontSize', FS, 'fontWeight', 'bold');
% $$$ text(timeSensor1(end)-18/1440, .75, sprintf('<\\gamma> = %1.2f',nanmean(JbVec./epsVec)),'fontSize',FS2,'fontWeight','bold')
% $$$ set(gca, 'tickdir', 'out')
% $$$ text(timeSensor1(1)+2/1440, .85, 'c','fontSize',FS2,'fontWeight', 'bold', 'BackgroundColor',[1 1 1]);
% $$$ set(gca, 'fontSize', FS)
% $$$ adjust_space


%% add Richardson number
POS = get(s1, 'pos');
ax = axes('pos', [0.08 POS(2) 0.12 POS(4)], 'yticklabel', [], 'xticklabel',[]);
plot(nanmean(Ri.Ri,2), totalDepth-Ri.zVec, 'k', 'linewidth', 2)
hold on
plot(Ri.zVec*0+.25, totalDepth-Ri.zVec, '--k')
set(ax, 'ydir', 'normal');
set(ax, 'xscale', 'log');
set(ax, 'tickdir', 'out')
ylim([0 max(totalDepth-zVecReg)])
ylab = ylabel('hab (m)', 'fontSize', FS, 'fontWeight', 'bold');
xlabel('Ri')
xlim([1e-2 1e4])
set(gca, 'xtick', [1e-2 1 1e2 1e4])
set(gca, 'fontSize', FS)
set(ylab, 'pos', [.0001, 60.4264, 1])

%outfile = sprintf('Reorder2D_Jb_20121009_%02dh25.png',myTimes(ii));
outfile = 'Jb-eps_yyyymmdd_.png';
set(gcf, 'renderer', 'painters')
print('-dpng', '-r300', outfile)


%outfile = sprintf('Jb-eps_%02d.eps',ii)
outfile = 'Jb-eps_toberename.eps';

print('-depsc2', outfile)


