

%% Plot section
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 5; % no. subplot row
dx = 0.015; % horiz. space between subplots
dy = 0.02; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.14; % very right of figure
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
s1 = subplot(511);
contourf(timeAve, zVec, tempAve, 100, 'linestyle', 'none')
set(gca, 'ydir', 'reverse')
set(gca, 'xticklabel', []);
cb = colorbar;
adjust_space
drawnow
hold on
for i = 1:length(timeAveHist)
    plot([timeAveHist(i) timeAveHist(i)],[zVec(1) zVec(end)], '--k')
end
hold off

cbPos = get(cb, 'pos');
cbPos(1) = cbPos(1)-dx;
cbPos(2) = cbPos(2)+dy;
cbPos(4) = cbPos(4)-2*dy;
cbPos(3) = cbPos(3)*.6;
set(cb, 'pos', cbPos);
ti = ylabel(cb,'T (^{\circ}C)', 'FontSize', 10);
drawnow
% $$$ ti_pos = get(ti, 'position');
% $$$ set(ti, 'Rotation',0.0);     
% $$$ clim = get(gca, 'clim');
% $$$ set(ti, 'position', [(clim(2)-clim(1))./2  2 ti_pos(3)]); 

%% s2
s2 = subplot(512);
%pcolor(timeAve, zVec, UAve); shading interp;
contourf(timeAve, zVec, UAve, 100, 'linestyle', 'none');
ylabel('Depth (m)')
set(gca, 'ydir', 'reverse')
set(gca, 'xticklabel', []);
%cb = colorbar
caxis([-.15 .15])
adjust_space
drawnow
hold on
for i = 1:length(timeAveHist)
    plot([timeAveHist(i) timeAveHist(i)],[zVec(1) zVec(end)], '--k')
end
hold off

%% s3
s3 = subplot(513);
%pcolor(timeAve, zVec, VAve); shading interp;
contourf(timeAve, zVec, VAve, 100, 'linestyle', 'none');
set(gca, 'ydir', 'reverse')
set(gca, 'xticklabel', []);
%set(cb, 'delete')
cb = colorbar;
caxis([-.15 .15])
adjust_space
drawnow
hold on
for i = 1:length(timeAveHist)
    plot([timeAveHist(i) timeAveHist(i)],[zVec(1) zVec(end)], '--k')
end
hold off

cbPos = get(cb, 'pos');
cbPos(1) = cbPos(1)-dx;
cbPos(2) = cbPos(2)+dy;
cbPos(4) = cbPos(4)*2-2*dy;
cbPos(3) = cbPos(3)*.6;
set(cb, 'pos', cbPos)
ti = ylabel(cb,'u,v (m s^{-1})', 'FontSize', 10);
drawnow

%% s4
s4 = subplot(514);
S2Ave(end-5:end,:) = NaN;
pcolor(timeAve, zVec, log10(S2Ave)); shading interp;
set(gca, 'ydir', 'reverse')
set(gca, 'xticklabel', []);
cb = colorbar;
caxis([-4.3 -3.5])
adjust_space
drawnow
hold on
for i = 1:length(timeAveHist)
    plot([timeAveHist(i) timeAveHist(i)],[zVec(1) zVec(end)], '--k')
end
hold off

cbPos = get(cb, 'pos');
cbPos(1) = cbPos(1)-dx;
cbPos(2) = cbPos(2)+dy;
cbPos(4) = cbPos(4)-2*dy;
cbPos(3) = cbPos(3)*.6;
set(cb, 'pos', cbPos)
ti = ylabel(cb,'log(S^2 / s^{-2})', 'FontSize', 10);
drawnow


%% s5
s5 = subplot(515);
pcolor(timeAve, eps.zVecReg, log10(epsAve)); shading interp
set(gca, 'ydir', 'reverse')
%xlabel('time to V_{max}')
xlabel('t (days)')
cb = colorbar;
adjust_space
drawnow
hold on
for i = 1:length(timeAveHist)
    plot([timeAveHist(i) timeAveHist(i)],[zVec(1) zVec(end)], '--k')
end
hold off

cbPos = get(cb, 'pos');
cbPos(1) = cbPos(1)-dx;
cbPos(2) = cbPos(2)+dy;
cbPos(4) = cbPos(4)-2*dy;
cbPos(3) = cbPos(3)*.6;
set(cb, 'pos', cbPos)
ti = ylabel(cb,'log(\epsilon / W kg^{-1})', 'FontSize', 10);
drawnow

set(gcf, 'renderer', 'painters')
print('-depsc2', 'tidal_average.eps')



figure(2)
pcolor(timeAve, zVec, log10(N2Ave)); shading flat;
set(gca, 'ydir', 'reverse')
%xlabel('time to V_{max}')
xlabel('t (days)')
cb = colorbar;
drawnow
hold on
for i = 1:length(timeAveHist)
    plot([timeAveHist(i) timeAveHist(i)],[zVec(1) zVec(end)], '--k')
end
hold off

cbPos = get(cb, 'pos');
cbPos(1) = cbPos(1)-dx;
cbPos(2) = cbPos(2)+dy;
cbPos(4) = cbPos(4)-2*dy;
cbPos(3) = cbPos(3)*.6;
set(cb, 'pos', cbPos)
ti = ylabel(cb,'log(N^{2} / s^{-2})', 'FontSize', 10);
drawnow