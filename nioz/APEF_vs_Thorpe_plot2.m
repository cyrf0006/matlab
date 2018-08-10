
%% Few plots
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 5; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.14; % very right of figure
tops = 0.03; % top of figure
bots = 0.07; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %
figure(1)
clf
set(gcf, 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 16 18])
s1 = subplot(5,1,[1 2]);
%imagesc(timeSensor1, zVecReg, Ritp1-1000)
imagesc(timeSensor1, zVecReg, Titp)
hold on
contour(timeSensor1, zVecReg, Titp, 10, 'color', 'k')
hold off
c = colorbar;
%ylabel(c, '\sigma_1 (kg m^{-3})', 'FontSize', 10); 
ylabel(c, 'T (^{\circ}C)', 'FontSize', 10); 
ylabel('Depth (m)')
%caxis([1032 1032.1]-1000)
datetick
set(gca, 'xticklabel', []);
xlim([timeSensor1(1) timeSensor1(end)])

adjust_space
pos1 = get(gca, 'pos');
adjust_space
pos2 = get(gca, 'pos');
pos2(4) = 2*pos2(4)+dy;
set(gca, 'pos', pos2)
drawnow

cbPos = get(c, 'pos');
cbPos(1) = cbPos(1)-.01;
cbPos(3) = cbPos(3)/2;
cbPos(2) = cbPos(2)+.05;
cbPos(4) = cbPos(4)-2*.05;
set(c, 'pos', cbPos)
drawnow

s2 = subplot(5,1,3);
semilogy(timeSensor1, JbVec, 'k') 
hold on
plot([timeSensor(1) timeSensor(end)], [1 1]*nanmean(JbVec), '--k')
semilogy(timeSensor1, epsVec, 'color', [1 1 1]*.5)
plot([timeSensor(1) timeSensor(end)], [1 1]*nanmean(epsVec), 'color',  [1 1 1]*.5, 'linestyle', '--')
%semilogy(timeSensor1, JbVec1, 'b') 
%semilogy(timeSensor1, epsVec1, '--b')
%semilogy(timeSensor1, JbVec2, 'k') 
%semilogy(timeSensor1, epsVec2, '--k')
datetick
set(gca, 'xticklabel', []);
xlim([timeSensor1(1) timeSensor1(end)])
%xlabel(sprintf(datestr(timeSensor(I(1)), 1)));
ylabel('J_b,\epsilon (m^2 s^{-3})')
set(gca, 'ytick', [1e-9 1e-8 1e-7 1e-6])
ylim([1e-9 5e-7])
adjust_space
drawnow


s3 = subplot(5,1,4);
%semilogy(timeSensor1, N2Vec, 'k') 
plot(timeSensor1, N2Vec, 'k') 
hold on
plot(timeSensor1, S2Vec, 'color',  [1 1 1]*.5) 

%lot([timeSensor(1) timeSensor(end)], [1 1]*nanmean(N2Vec), '--k')
%emilogy(timeSensor1, S2Vec, 'color', [1 1 1]*.5)
%lot([timeSensor(1) timeSensor(end)], [1 1]*nanmean(S2Vec), 'color',  [1 1 1]*.5, 'linestyle', '--')
datetick
set(gca, 'xticklabel', []);
xlim([timeSensor1(1) timeSensor1(end)])
%xlabel(sprintf(datestr(timeSensor(I(1)), 1)));
ylabel('N^2,S^2 (s^{-2})')
set(gca, 'ytick', [1e-9 1e-8 1e-7 1e-6])
%lim([1e-9 5e-7])
adjust_space
drawnow



s4 = subplot(5,1,5);
plot(timeSensor1,JbVec./epsVec, 'k')
%datetick('x', 7)
hold on
plot(timeSensor1,N2Vec./S2Vec, 'color', [1 1 1]*.5)
%plot(timeSensor1,JbVec1./epsVec1, 'b')
%plot(timeSensor1,JbVec2./epsVec2, 'k')
plot(timeSensor1,JbVec*0+.2, 'color', [1 1 1]*.0,'linestyle', '--')
plot(timeSensor1,JbVec*0+.25, 'color', [1 1 1]*.5,'linestyle', '--')
ylim([0 1])
datetick('x', 15)
xlim([timeSensor1(1) timeSensor1(end)])
xlabel(sprintf(datestr(timeSensor1(1), 1)));
%xlabel('2012-10-11')
%xlabel('October 2012')
ylabel('\Gamma');
text(timeSensor1(end)-10/1440, .75, sprintf('<\\Gamma> = %1.2f',nanmean(JbVec./epsVec)))
%text(timeSensor1(end)-10/1440, .60, sprintf('Gam_2 = %1.2f',nanmean(JbVec1./epsVec1)))
%text(timeSensor1(end)-10/1440, 1.25, sprintf('Gam_3 = %1.2f',nanmean(JbVec2./epsVec2)))
adjust_space

%outfile = sprintf('Reorder2D_Jb_20121009_%02dh25.png',myTimes(ii));
outfile = 'Jb-eps_yyyymmdd_.png';
set(gcf, 'renderer', 'painters')
print('-dpng', '-r300', outfile)

outfile = sprintf('Jb-eps_%02d.eps',ii)
print('-depsc2', outfile)

