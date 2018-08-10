
%% Few plots
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 2; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.12; % very left of figure
rigs = 0.14; % very right of figure
tops = 0.03; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

x = [-.125 .125 .125 -.125 -.125];
x2 = timeAve;

figure(1)
clf
set(gcf, 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 16 12])

%% S1
s1 = subplot(2,1,1);
YLIM = [1e-8 1e-6];
patch(x, [YLIM(1) YLIM(1) YLIM(2) YLIM(2) YLIM(1)], [1 1 1]*.8, 'linestyle', 'none');
set(gca, 'yscale', 'log')
hold on
%YLIM = get(gca, 'ylim')
%patch(x, [YLIM(1) YLIM(1) YLIM(2) YLIM(2) YLIM(1)], [1 1 1]*.8, 'linestyle', 'none');

yy1 = eps2p5;
yy2 = eps97p5;
patch([x2 fliplr(x2) x2(1)], [yy1 fliplr(yy2) yy1(1)], [.5 .5 .5], 'linestyle', 'none')
semilogy(timeAve, epsAve, 'color', [.5 .5 .5], 'linewidth', 2)
semilogy(timeAve, eps2p5, 'color', [.5 .5 .5], 'linewidth', 2)
semilogy(timeAve, eps97p5, 'color', [.5 .5 .5], 'linewidth', 2)

yy1 = Jb2p5;
yy2 = Jb97p5;
patch([x2 fliplr(x2) x2(1)], [yy1 fliplr(yy2) yy1(1)], [0 0 0], 'linestyle', 'none')
semilogy(timeAve, JbAve, 'color', [.5 .5 .5]*0, 'linewidth', 2)
semilogy(timeAve, Jb2p5, 'color', [.5 .5 .5]*0, 'linewidth', 2)
semilogy(timeAve, Jb97p5, 'color', [.5 .5 .5]*0, 'linewidth', 2)
ylim(YLIM);
set(gca, 'xticklabel', '')
set(gca, 'box', 'on')
set(gca, 'tickdir', 'out')
set(gca, 'xtick', [-.5:.2:.5])
set(gca, 'xgrid', 'on')
ylabel('\epsilon, J_b (W kg^{-1})', 'fontSize', FS)
hold off
text(.45, 1.5e-8, 'a', 'fontSize', FS2, 'fontWeight', 'bold')
text(.52, 1e-7, 'J_b', 'fontSize', FS2, 'fontWeight', 'bold')
text(.52, 8e-7, '\epsilon', 'fontSize', FS2, 'fontWeight', 'bold')
set(gca, 'fontSize', FS)
adjust_space
box on

%% S2
s2 = subplot(2, 1,2);
YLIM = sqrt([1e-6 1e-4]);
patch(x, [YLIM(1) YLIM(1) YLIM(2) YLIM(2) YLIM(1)], [1 1 1]*.8, 'linestyle', 'none');
set(gca, 'yscale', 'log')

hold on
yy1 =sqrt(N22p5);
yy2 = sqrt(N297p5);
patch([x2 fliplr(x2) x2(1)], [yy1 fliplr(yy2) yy1(1)], [0 0 0], ...
      'linestyle', 'none')

ylim(YLIM);
set(gca, 'box', 'on')
set(gca, 'tickdir', 'out')
set(gca, 'xtick', [-.5:.2:.5])
set(gca, 'xgrid', 'on')
ylabel('N (s^{-1})', 'fontSize', FS)
semilogy(timeAve, sqrt(N2Ave), 'k', 'linewidth', 2)
semilogy(timeAve, sqrt(N22p5), 'k', 'linewidth', 2)
semilogy(timeAve, sqrt(N297p5), 'k', 'linewidth', 2)
hold off
text(.45, 1.2e-3, 'b', 'fontSize', FS2, 'fontWeight', 'bold')
set(gca, 'fontSize', FS)
xlabel('\phi (day)', 'fontSize', FS)
adjust_space
box on



I = find(timeAve>=-.125 & timeAve<= .125);
II = find(timeAve<-.125 | timeAve> .125);

[nanmean(epsAve(I)) nanmean(epsAve(II)) nanmean(epsAve(II))/nanmean(epsAve(I))]
[nanmean(JbAve(I)) nanmean(JbAve(II)) nanmean(JbAve(I))/nanmean(JbAve(II))]
[nanmean(sqrt(N2Ave(I))) nanmean(sqrt(N2Ave(II))) nanmean(sqrt(N2Ave(II)))/nanmean(sqrt(N2Ave(I)))]

set(gcf, 'renderer', 'painters')
outfile = 'epsTidal.eps';
print('-depsc2', outfile)
