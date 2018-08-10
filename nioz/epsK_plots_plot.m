figure(2)
clf
set(gcf, 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 16 6])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.12; % very left of figure
rigs = 0.05; % very right of figure
tops = 0.03; % top of figure
bots = 0.2; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

%% S1
x = [-.25 .2 .2 -.25 -.25];
YLIM = [1500 3250];
patch(x, [YLIM(1) YLIM(1) YLIM(2) YLIM(2) YLIM(1)], [1 1 1]*.8, 'linestyle', 'none');
hold on
yy1 = 2*pi./sqrt(N22p5);
yy2 = 2*pi./sqrt(N297p5);
patch([x2 fliplr(x2) x2(1)], [yy1 fliplr(yy2) yy1(1)], [.5 .5 .5], 'linestyle', 'none')
yy1 = TN2p5;
yy2 = TN97p5;
patch([x2 fliplr(x2) x2(1)], [yy1 fliplr(yy2) yy1(1)], [0 0 0], 'linestyle', 'none')
set(gca, 'box', 'on')
set(gca, 'tickdir', 'out')
set(gca, 'xtick', [-.5:.2:.5])
set(gca, 'xgrid', 'on')
set(gca, 'ygrid', 'on')
ylabel('\tau (s)', 'fontSize', FS)

yy1 = 2*pi./sqrt(S2bg2p5);
yy2 = 2*pi./sqrt(S2bg97p5);
patch([x2 fliplr(x2) x2(1)], [yy1 fliplr(yy2) yy1(1)], [.5 .5 1], 'linestyle', 'none')

plot(x2, 2*pi./sqrt(S2bgAve), 'color', [.5 .5 1])

hold off
text(.05, 3000, '\tau_N', 'fontSize', FS2, 'fontWeight', 'bold')
text(.3, 2300, '\tau_C', 'fontSize', FS2, 'fontWeight', 'bold')
%text(.4, 800, '\tau_b', 'fontSize', FS2, 'fontWeight', 'bold')
set(gca, 'fontSize', FS)
adjust_space
xlabel('\phi (day)', 'fontSize', FS)
ylim([YLIM(1)-1 YLIM(end)+1])
box on




set(gcf, 'renderer', 'painters')
outfile = 'tauTidal.eps';
print('-depsc2', outfile)

% $$$ outfile = 'tauTidal_epsutils.eps';
% $$$ epswrite(outfile, 'EmbedFonts', 'all')

