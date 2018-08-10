FS = 12;
FS2 = 15;
close all
%% Plot section
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 6; % no. subplot row
dx = 0.015; % horiz. space between subplots
dy = 0.02; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.15; % very right of figure
tops = 0.03; % top of figure
bots = 0.07; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

%% Colormap
load BR_symetric
%colormap(BR_symetric)
cm = colormap(jet);
% $$$ cm = cm*1.3;
% $$$ I = find(cm > 1);
% $$$ cm(I) = 1;
% or
cm = cm*.75;
cm = colormap(cm);


figure(1)
clf
set(gcf, 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 16 20])
%%s1
s1 = subplot(611);
contourf(timeAve, totalDepth-zVec, tempAve, 100, 'linestyle', 'none')
set(gca, 'ydir', 'reverse')
set(gca, 'xticklabel', []);
set(gca, 'tickdir', 'out');
set(gca, 'ydir', 'normal');
cb = colorbar;
adjust_space
drawnow
hold on
for i = 1:length(timeAveHist)
    plot([timeAveHist(i) timeAveHist(i)],[totalDepth-zVec(1) totalDepth-zVec(end)], '--k')
end
hold off
ylabel('hab (m)', 'fontSize', FS, 'fontWeight', 'bold')
text(.45, 24, 'a', 'fontSize', FS2, 'fontWeight', 'bold','color', [1 1 1])
ylim([0 120])

cbPos = get(cb, 'pos');
cbPos(1) = cbPos(1)-dx;
cbPos(2) = cbPos(2)+dy;
cbPos(4) = cbPos(4)-2*dy;
cbPos(3) = cbPos(3)*.6;
set(cb, 'pos', cbPos);
set(gca, 'fontSize', FS)
ti = ylabel(cb,'\theta_0(^{\circ}C)', 'FontSize', FS, 'fontWeight', 'bold')
drawnow
% $$$ ti_pos = get(ti, 'position');
% $$$ set(ti, 'Rotation',0.0);     
% $$$ clim = get(gca, 'clim');
% $$$ set(ti, 'position', [(clim(2)-clim(1))./2  2 ti_pos(3)]); 

%% s2
s2 = subplot(612);
%pcolor(timeAve, totalDepth-zVec, UAve); shading interp;
contourf(timeAve, totalDepth-zVec, UAve, 100, 'linestyle', 'none');
ylabel('hab (m)', 'fontSize', FS, 'fontWeight', 'bold')
set(gca, 'ydir', 'reverse')
set(gca, 'xticklabel', []);
set(gca, 'tickdir', 'out');
set(gca, 'ydir', 'normal');
set(gca, 'fontSize', FS)
cb = colorbar;
caxis([-.15 .15])
adjust_space
drawnow
hold on
for i = 1:length(timeAveHist)
    plot([timeAveHist(i) timeAveHist(i)],[totalDepth-zVec(1) totalDepth-zVec(end)], '--k')
end
hold off
text(.45, 24, 'b', 'fontSize', FS2, 'fontWeight', 'bold','color', [1 1 1])
ylim([0 120])
cbPos = get(cb, 'pos');
cbPos(1) = cbPos(1)-dx;
cbPos(2) = cbPos(2)+dy;
cbPos(4) = cbPos(4)-2*dy;
cbPos(3) = cbPos(3)*.6;
set(cb, 'pos', cbPos)
ti = ylabel(cb,'u (m s^{-1})', 'FontSize', FS, 'fontWeight', 'bold');
set(gca, 'fontSize', FS)
drawnow

%% s3
s3 = subplot(613);
%pcolor(timeAve, totalDepth-zVec, VAve); shading interp;
contourf(timeAve, totalDepth-zVec, VAve, 100, 'linestyle', 'none');
ylabel('hab (m)', 'fontSize', FS, 'fontWeight', 'bold')
set(gca, 'ydir', 'reverse')
set(gca, 'xticklabel', []);
set(gca, 'tickdir', 'out');
set(gca, 'ydir', 'normal');
%set(cb, 'delete')
cb = colorbar;
caxis([-.15 .15])
adjust_space
drawnow
hold on
for i = 1:length(timeAveHist)
    plot([timeAveHist(i) timeAveHist(i)],[totalDepth-zVec(1) totalDepth-zVec(end)], '--k')
end
hold off
text(.45, 24, 'c', 'fontSize', FS2, 'fontWeight', 'bold','color', [1 1 1])
ylim([0 120])

cbPos = get(cb, 'pos');
cbPos(1) = cbPos(1)-dx;
cbPos(2) = cbPos(2)+dy;
cbPos(4) = cbPos(4)-2*dy;
cbPos(3) = cbPos(3)*.6;
set(cb, 'pos', cbPos)
ti = ylabel(cb,'v (m s^{-1})', 'FontSize', FS, 'fontWeight', 'bold');
set(gca, 'fontSize', FS)

% $$$ cbPos = get(cb, 'pos');
% $$$ cbPos(1) = cbPos(1)-dx;
% $$$ cbPos(2) = cbPos(2)+dy;
% $$$ cbPos(4) = cbPos(4)*2-2*dy;
% $$$ cbPos(3) = cbPos(3)*.6;
% $$$ set(cb, 'pos', cbPos)
% $$$ set(cb, 'ytick', [-.15:.05:.15])
% $$$ set(cb, 'yticklabel', [-.15:.05:.15])
% $$$ ti = ylabel(cb,'u,v (m s^{-1})', 'FontSize', FS);
% $$$ set(gca, 'fontSize', FS)

drawnow

%% s4
s4 = subplot(614);
S2Ave(end-5:end,:) = NaN;
%contourf(timeAve, totalDepth-Ri.zVec, log10(S2Ave), 100, 'linestyle', 'none');
pcolor(timeAve, totalDepth-Ri.zVec, log10(S2Ave)); shading interp;
ylabel('hab (m)', 'fontSize', FS, 'fontWeight', 'bold')
set(gca, 'ydir', 'reverse')
set(gca, 'xticklabel', []);
set(gca, 'tickdir', 'out');
set(gca, 'ydir', 'normal');
cb = colorbar;
%caxis([-5.4 -4.6])
%caxis([-5.2 -4.2])
adjust_space
drawnow
hold on
for i = 1:length(timeAveHist)
    plot([timeAveHist(i) timeAveHist(i)],[totalDepth-zVec(1) totalDepth-zVec(end)], '--k')
end
hold off
text(.45, 24, 'd', 'fontSize', FS2, 'fontWeight', 'bold','color', [1 1 1])
ylim([0 120])

cbPos = get(cb, 'pos');
cbPos(1) = cbPos(1)-dx;
cbPos(2) = cbPos(2)+dy;
cbPos(4) = cbPos(4)-2*dy;
cbPos(3) = cbPos(3)*.6;
set(cb, 'pos', cbPos)
ti = ylabel(cb,'log(S^2 / s^{-2})', 'FontSize', FS, 'fontWeight', 'bold');
set(gca, 'fontSize', FS)
drawnow

%% s5
s5 = subplot(615);
S2Ave(end-5:end,:) = NaN;
%contourf(timeAve, totalDepth-Ri.zVec, log10(N2Ave), 100, 'linestyle', 'none');
pcolor(timeAve, totalDepth-Ri.zVec, log10(N2Ave)); shading interp;
ylabel('hab (m)', 'fontSize', FS, 'fontWeight', 'bold')
set(gca, 'ydir', 'reverse')
set(gca, 'xticklabel', []);
set(gca, 'tickdir', 'out');
set(gca, 'ydir', 'normal');
cb = colorbar;
%caxis([-5 -4.2])
%caxis([-5 -4.6])
%caxis([-5.2 -4.2])
adjust_space
drawnow
hold on
for i = 1:length(timeAveHist)
    plot([timeAveHist(i) timeAveHist(i)],[totalDepth-zVec(1) totalDepth-zVec(end)], '--k')
end
hold off
text(.45, 24, 'e', 'fontSize', FS2, 'fontWeight', 'bold','color', [1 1 1])
ylim([0 120])

% $$$ % For double colorbar
% $$$ cbPos = get(cb, 'pos');
% $$$ cbPos(1) = cbPos(1)-dx;
% $$$ cbPos(2) = cbPos(2)+dy;
% $$$ cbPos(4) = cbPos(4)*2-2*dy;
% $$$ cbPos(3) = cbPos(3)*.6;

cbPos = get(cb, 'pos');
cbPos(1) = cbPos(1)-dx;
cbPos(2) = cbPos(2)+dy;
cbPos(4) = cbPos(4)-2*dy;
cbPos(3) = cbPos(3)*.6;
set(cb, 'pos', cbPos)
ti = ylabel(cb,'log(N^2 / s^{-2})', 'FontSize', FS, 'fontWeight', 'bold');
set(gca, 'fontSize', FS)
drawnow


%% s6
s6 = subplot(616);
%pcolor(timeAve, eps.totalDepth-zVecReg, log10(epsAve)); shading interp
contourf(timeAve, totalDepth-eps.zVecReg, log10(epsAve), 100, 'linestyle', 'none');
ylabel('hab (m)', 'fontSize', FS, 'fontWeight', 'bold')
set(gca, 'ydir', 'reverse')
set(gca, 'tickdir', 'out');
set(gca, 'ydir', 'normal');
set(gca, 'xtick', [-.5:.2:.5])
%xlabel('time to V_{max}')
xlabel('\phi (days)', 'FontSize', FS, 'fontWeight', 'bold')
cb = colorbar;
adjust_space
drawnow
hold on
for i = 1:length(timeAveHist)
    plot([timeAveHist(i) timeAveHist(i)],[totalDepth-zVec(1) totalDepth-zVec(end)], '--k')
end
hold off
text(.45, 24, 'f', 'fontSize', FS2, 'fontWeight', 'bold','color', [1 1 1])
ylim([0 120])

cbPos = get(cb, 'pos');
cbPos(1) = cbPos(1)-dx;
cbPos(2) = cbPos(2)+dy;
cbPos(4) = cbPos(4)-2*dy;
cbPos(3) = cbPos(3)*.6;
set(cb, 'pos', cbPos)
ti = ylabel(cb,'log(\epsilon / W kg^{-1})', 'FontSize', FS, 'fontWeight', 'bold');
set(gca, 'fontSize', FS)
drawnow

set(gcf, 'renderer', 'painters')
print('-depsc2', 'tidal_average.eps')


