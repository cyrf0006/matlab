xbot = [timeVec; timeVec(end); timeVec(1); timeVec(1)]; 
ybot = [botFilt; zf; zf; botFilt(1)];
FS1 = 14;
FS2 = 16;
FS0 = 12;


figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 25 20])
pcolor(timeVec, zVec, TMat)
shading interp
hold on
contour(timeVec, zVec, RMat, [20:.5:30], 'lineColor', 'k')
keyboard
patch(xbot, ybot, [1 .9333 .6667]/1.5)
drawnow
caxis([13 20])
set(gca, 'ydir', 'reverse')
set(gca, 'xdir', 'reverse')
datetick('x', 7)
colorbar
xlim([t0 tf])
xlabel('Nov. 2014', 'fontSize', FS1, 'fontWeight', 'bold')
ylabel('Depth (m)', 'fontSize', FS1, 'fontWeight', 'bold')
title('Temperature (^{\circ}C)', 'fontSize', FS2, 'fontWeight', 'bold')
set(gca, 'fontSize', FS0, 'fontWeight', 'bold')
print('-dpng', '-r300', 'M241_temp.png')

figure(2)
clf
pcolor(timeVec, zVec, SMat)
shading interp
hold on
contour(timeVec, zVec, RMat, [20:.5:30], 'lineColor', 'k')
caxis([34 39])
set(gca, 'ydir', 'reverse')
set(gca, 'xdir', 'reverse')
datetick
colorbar


figure(3)
clf
pcolor(timeVec, zVec, chlMat)
shading interp
hold on
contour(timeVec, zVec, RMat, [20:.5:30], 'lineColor', 'k')
%caxis([34 39])
set(gca, 'ydir', 'reverse')
set(gca, 'xdir', 'reverse')
colorbar
datetick
caxis([0 1])

figure(4)
clf
pcolor(timeVec, zVec, cdMat)
shading interp
hold on
contour(timeVec, zVec, RMat, [20:.5:30], 'lineColor', 'k')
caxis([1.5 4])
set(gca, 'ydir', 'reverse')
set(gca, 'xdir', 'reverse')
colorbar
datetick
%caxis([0 1])

figure(5)
clf
pcolor(timeVec, zVec, log10(bbMat))
shading interp
hold on
contour(timeVec, zVec, RMat, [20:.5:30], 'lineColor', 'k')
%caxis([34 39])
set(gca, 'ydir', 'reverse')
set(gca, 'xdir', 'reverse')
colorbar
datetick
caxis([-5 -2])