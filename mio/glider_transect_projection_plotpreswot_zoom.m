%% THIS SCRIPT GENERATES THE FIGURE FOUND IN MANUSCRIPT ON PAH
% -> Run glider_proj_preswot first!

%% Manual removal of problematic counts in NAP
% $$$ NAP(:,36) = NaN;
% $$$ NAP(75:300,92) = NaN;
% $$$ NAP(1:40,72) = NaN;
% $$$ NAP(1:100,74) = NaN;
zMax = 300;

xtra_offset = 0.03;
cbar_offset = 0.02;
FONTSIZE = 10;
vredux = .8;
hredux = .6;

textx = 5;
texty = 275;
texty3 = 20;

FS1=12; % label
%XLIM = [0 m_lldist([origin(2) target(2)], [origin(1) target(1)],1)]; %<-- aimed transect
XLIM = [0 105];
XTICKS = 0:20:XLIM(end);
v1 = [27:.1:29];
v2 = [27:.2:29];

%% Save info:
timeProj = timeVec(theIndex);
latProj = latVec(theIndex);
lonProj = lonVec(theIndex);
save projection_nexos02.mat timeProj latProj lonProj


x1 = [1 20]; %km
x2 = [45 55];
x3 = [70 80];
xBox1 = [x1(1) x1(2) x1(2) x1(1) x1(1)];
xBox2 = [x2(1) x2(2) x2(2) x2(1) x2(1)];
xBox3 = [x3(1) x3(2) x3(2) x3(1) x3(1)];
yBox = [zMin zMin zMax zMax zMin];

figure(1) % Temp
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10 5.5])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.04; % horiz. space between subplots
dy = 0.1; % vert. space between subplots
lefs = 0.15; % very left of figure
rigs = 0.12; % very right of figure
tops = 0.02; % top of figure
bots = 0.16; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;

count_col = 1;
count_row = 1;
% *************************************************************** %

% S1
pcolor(xVecCTD,zVec, CT(:,theIndexCTD)); shading interp
hold on
[C,h] = contour(xVecCTD, zVec, sig0(:,theIndexCTD), v1, 'color', [.5 .5 .5]);
clabel(C,h,v2, 'color', [.5 .5 .5])
map = flipud(brewermap(14, 'RdBu'));
colormap(map);
caxis([13 17])
%plot(xBox3, yBox, '--','color', [1 1 1]*0, 'lineWidth', 2)
ylabel('Depth (m)', 'FontSize', FS1, 'fontWeight', 'bold')
set(gca, 'ygrid', 'on');
set(gca, 'ydir', 'reverse');
set(gca, 'xtick', XTICKS)
set(gca, 'tickdir', 'out');
xlim(XLIM)
ylim([0 zMax])
text(XLIM(1)+2,zMax-.1*zMax, 'T (^{\circ}C)  (A)', 'horizontalAlignment', 'left', 'fontSize', FS1, 'fontWeight', 'bold')
cb = colorbar;
adjust_space
pause(1)

cpos = get(cb, 'pos');
cpos(1) = cpos(1) - cbar_offset;
cpos(2) = cpos(2)+(cpos(4)-cpos(4)*vredux)/2;
cpos(4) = cpos(4)*vredux;
cpos(3) = cpos(3)*hredux;
set(cb, 'pos', cpos)
set(gcf, 'renderer', 'painters')
print(gcf, '-dpng', '-r300', 'M327_T_noplot.png')  


%% --------------------------------------------------------------------- %

figure(2) % Sal
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10 5.5])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.04; % horiz. space between subplots
dy = 0.1; % vert. space between subplots
lefs = 0.11; % very left of figure
rigs = 0.15; % very right of figure
tops = 0.02; % top of figure
bots = 0.16; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

pcolor(xVecCTD,zVec, SA(:,theIndexCTD)); shading interp
hold on
[C,h] = contour(xVecCTD, zVec, sig0(:,theIndexCTD), v1, 'color', [.5 .5 .5]);
clabel(C,h,v2, 'color', [.5 .5 .5])
map = flipud(brewermap(14, 'RdBu'));
colormap(map);
caxis([37.5 39])
%plot(xBox3, yBox, '--','color', [1 1 1]*0, 'lineWidth', 2)
set(gca, 'yticklabel', []);
set(gca, 'ygrid', 'on');
set(gca, 'ydir', 'reverse');
set(gca, 'xtick', XTICKS)
set(gca, 'tickdir', 'out');
xlim(XLIM)
ylim([0 zMax])
%xlabel('Along-transect distace (km)', 'FontSize', 10, 'fontWeight', 'bold')
cb = colorbar;
%ti = ylabel(cb,'S_A (g Kg^{-1})', 'FontSize', 10, 'fontweight', 'bold');
text(XLIM(1)+2,zMax-.1*zMax, 'S_A (g Kg^{-1})  (B)', 'horizontalAlignment', 'left', 'fontSize', FS1, 'fontWeight', 'bold', 'color',[1 1 1]*0)
adjust_space
pause(1)

cpos = get(cb, 'pos');
cpos(1) = cpos(1) - cbar_offset;
cpos(2) = cpos(2)+(cpos(4)-cpos(4)*vredux)/2;
cpos(4) = cpos(4)*vredux;
cpos(3) = cpos(3)*hredux;
set(cb, 'pos', cpos)
set(gcf, 'renderer', 'painters')
print(gcf, '-dpng', '-r300', 'M327_S_noplot.png')  



%% --------------------------------------------------------------------- %

figure(3) % CHL
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10 5.5])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.04; % horiz. space between subplots
dy = 0.1; % vert. space between subplots
lefs = 0.11; % very left of figure
rigs = 0.15; % very right of figure
tops = 0.02; % top of figure
bots = 0.16; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

pcolor(xVec, zVec, CHL(:,theIndex)); shading flat
hold on
contour(xVecCTD, zVec, sig0(:,theIndexCTD), v1, 'color', [.5 .5 .5])
load('PuBuGn_modif.mat');
colormap(PuBuGn_modif);
caxis([0 0.8])
set(gca, 'yticklabel', []);
set(gca, 'ygrid', 'on');
set(gca, 'ydir', 'reverse');
set(gca, 'xtick', XTICKS)
set(gca, 'tickdir', 'out');
xlim(XLIM)
ylim([0 zMax])
cb = colorbar;
text(XLIM(1)+2,zMax-.1*zMax, '[Chl-a] (\mug L^{-1})  (D)', 'horizontalAlignment', 'left', 'fontSize', FS1, 'fontWeight', 'bold')
adjust_space
pause(1)

cpos = get(cb, 'pos');
cpos(1) = cpos(1) - cbar_offset;
cpos(2) = cpos(2)+(cpos(4)-cpos(4)*vredux)/2;
cpos(4) = cpos(4)*vredux;
cpos(3) = cpos(3)*hredux;
set(cb, 'pos', cpos)
set(gcf, 'renderer', 'painters')
print(gcf, '-dpng', '-r300', 'M327_CHL_noplot.png')  



%% --------------------------------------------------------------------- %
figure(4) % NAP
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10 5.5])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.04; % horiz. space between subplots
dy = 0.1; % vert. space between subplots
lefs = 0.15; % very left of figure
rigs = 0.12; % very right of figure
tops = 0.02; % top of figure
bots = 0.16; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

pcolor(xVec, zVec, TRYru(:,theIndex)); shading interp
hold on
contour(xVecCTD, zVec, sig0(:,theIndexCTD), v1, 'color', [.5 .5 .5])
% $$$ map = brewermap(14, 'PuBuGn'); 
% $$$ colormap(map);
load('PuBuGn_modif.mat');
colormap(PuBuGn_modif);
caxis([0.06, .15])
%caxis([0, .3]) %<-- for concentration
ylabel('Depth (m)', 'FontSize', FS1, 'fontWeight', 'bold')
set(gca, 'ygrid', 'on');
set(gca, 'ydir', 'reverse');
set(gca, 'xtick', XTICKS)
set(gca, 'tickdir', 'out');
xlim(XLIM)
ylim([0 zMax])
cb = colorbar;
text(XLIM(1)+2,zMax-.1*zMax, '[Naph-like] (RU)  (e)', 'horizontalAlignment', 'left', 'fontSize', FS1, 'fontWeight', 'bold')

adjust_space
pause(1)

cpos = get(cb, 'pos');
cpos(1) = cpos(1) - cbar_offset;
cpos(2) = cpos(2)+(cpos(4)-cpos(4)*vredux)/2;
cpos(4) = cpos(4)*vredux;
cpos(3) = cpos(3)*hredux;
set(cb, 'pos', cpos)
set(gcf, 'renderer', 'painters')
print(gcf, '-dpng', '-r300', 'M327_TRY_noplot.png')  


%% --------------------------------------------------------------------- %


figure(5) % CDOM
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10 5.5])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.04; % horiz. space between subplots
dy = 0.1; % vert. space between subplots
lefs = 0.15; % very left of figure
rigs = 0.12; % very right of figure
tops = 0.02; % top of figure
bots = 0.2; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

% S1
pcolor(xVec,zVec, CDOM(:,theIndex)); shading interp
hold on
contour(xVecCTD, zVec, sig0(:,theIndexCTD), v1, 'color', [.5 .5 .5])
map = brewermap(16, 'YlOrBr');        
colormap(map)
caxis([.9 1.5])
ylabel('Depth (m)', 'FontSize', FS1, 'fontWeight', 'bold')
set(gca, 'ygrid', 'on');
set(gca, 'ydir', 'reverse');
set(gca, 'xtick', XTICKS)
set(gca, 'tickdir', 'out');
xlim(XLIM)
ylim([0 zMax])
xlabel('Along-transect distace (km)', 'FontSize', 10, 'fontWeight', 'bold')
cb = colorbar;
text(XLIM(1)+2,zMax-.1*zMax, '[CDOM] (\mug L^{-1})  (G)', 'horizontalAlignment', 'left', 'fontSize', FS1, 'fontWeight', 'bold')

adjust_space
pause(1)

cpos = get(cb, 'pos');
cpos(1) = cpos(1) - cbar_offset;
cpos(2) = cpos(2)+(cpos(4)-cpos(4)*vredux)/2;
cpos(4) = cpos(4)*vredux;
cpos(3) = cpos(3)*hredux;
set(cb, 'pos', cpos)
set(gcf, 'renderer', 'painters')
pause(.5)
print(gcf, '-dpng', '-r300', 'M327_CDOM_noplot.png')  



%% --------------------------------------------------------------------- %
figure(6) % PHE
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10 5.5])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.04; % horiz. space between subplots
dy = 0.1; % vert. space between subplots
lefs = 0.11; % very left of figure
rigs = 0.15; % very right of figure
tops = 0.02; % top of figure
bots = 0.16; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

%% ATTENTION PHE TRES SENSIBLE AU BLANK...
pcolor(xVec, zVec, PHEru(:,theIndex)); shading interp
hold on
contour(xVecCTD, zVec, sig0(:,theIndexCTD), v1, 'color', [.5 .5 .5])
load('PuBuGn_modif.mat');
colormap(PuBuGn_modif);
caxis([0.03, .05])
%caxis([15, 22]) %<--- for concentrations
set(gca, 'ygrid', 'on');
set(gca, 'ydir', 'reverse');
set(gca, 'xtick', XTICKS)
set(gca, 'tickdir', 'out');
set(gca, 'yticklabel', [])
xlim(XLIM)
ylim([0 zMax])
cb = colorbar;
text(XLIM(1)+2,zMax-.1*zMax, '[Phe-like] (RU)  (f)', 'horizontalAlignment', 'left', 'fontSize', FS1, 'fontWeight', 'bold')

adjust_space
pause(1)

cpos = get(cb, 'pos');
cpos(1) = cpos(1) - cbar_offset;
cpos(2) = cpos(2)+(cpos(4)-cpos(4)*vredux)/2;
cpos(4) = cpos(4)*vredux;
cpos(3) = cpos(3)*hredux;
set(cb, 'pos', cpos)
set(gcf, 'renderer', 'painters')
print(gcf, '-dpng', '-r300', 'M327_PHE_noplot.png')  



%% --------------------------------------------------------------------- %


figure(7) % Turbidity
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10 5.5])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.04; % horiz. space between subplots
dy = 0.1; % vert. space between subplots
lefs = 0.11; % very left of figure
rigs = 0.15; % very right of figure
tops = 0.02; % top of figure
bots = 0.2; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

pcolor(xVec, zVec, log10(BB(:,theIndex))); shading interp
hold on
contour(xVecCTD, zVec, sig0(:,theIndexCTD), v1, 'color', [.5 .5 .5])
map = brewermap(16, 'YlOrBr');        
colormap(map)
caxis([-4 -3.5])
set(gca, 'yticklabel', []);
set(gca, 'ygrid', 'on');
set(gca, 'ydir', 'reverse');
set(gca, 'xtick', XTICKS)
set(gca, 'tickdir', 'out');
xlim(XLIM)
ylim([0 zMax])
xlabel('Along-transect distace (km)', 'FontSize', 10, 'fontWeight', 'bold')
cb = colorbar;
%ti = ylabel(cb,'log_{10}(Backscatter)', 'FontSize', 10, 'fontweight', 'bold');
text(XLIM(1)+2,zMax-.1*zMax, 'log_{10}(BB700)  (H)', 'horizontalAlignment', 'left', 'fontSize', FS1, 'fontWeight', 'bold')

adjust_space
pause(1)

cpos = get(cb, 'pos');
cpos(1) = cpos(1) - cbar_offset;
cpos(2) = cpos(2)+(cpos(4)-cpos(4)*vredux)/2;
cpos(4) = cpos(4)*vredux;
cpos(3) = cpos(3)*hredux;
set(cb, 'pos', cpos)
set(gcf, 'renderer', 'painters')
print(gcf, '-dpng', '-r300', 'M327_BB_noplot.png')  





figure(8) % N2
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10 5.5])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.04; % horiz. space between subplots
dy = 0.1; % vert. space between subplots
lefs = 0.15; % very left of figure
rigs = 0.12; % very right of figure
tops = 0.02; % top of figure
bots = 0.16; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

pcolor(xVecCTD, zVec, log10(N2(:,theIndexCTD))); shading interp
hold on
contour(xVecCTD, zVec, sig0(:,theIndexCTD), v1, 'color', [.5 .5 .5])
map = flipud(brewermap(14, 'RdBu'));
colormap(map)
caxis([-4.7 -3.2])
ylabel('Depth (m)', 'FontSize', FS1, 'fontWeight', 'bold')
set(gca, 'ygrid', 'on');
set(gca, 'ydir', 'reverse');
set(gca, 'xtick', XTICKS)
set(gca, 'tickdir', 'out');
xlim(XLIM)
ylim([0 zMax])
%xlabel('Along-transect distace (km)', 'FontSize', 10, 'fontWeight', 'bold')
cb = colorbar;
%ti = ylabel(cb,'log_{10}(N^2)', 'FontSize', 10, 'fontweight', 'bold');
text(XLIM(1)+2,zMax-.1*zMax, 'log_{10}(N^2 (s^{-2}))  (C)', 'horizontalAlignment', 'left', 'fontSize', FS1, 'fontWeight', 'bold')

adjust_space
pause(1)

cpos = get(cb, 'pos');
cpos(1) = cpos(1) - cbar_offset;
cpos(2) = cpos(2)+(cpos(4)-cpos(4)*vredux)/2;
cpos(4) = cpos(4)*vredux;
cpos(3) = cpos(3)*hredux;
set(cb, 'pos', cpos)
set(gcf, 'renderer', 'painters')
print(gcf, '-dpng', '-r300', 'M327_N2_noplot.png')  




figure(9) % O2
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10 5.5])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.04; % horiz. space between subplots
dy = 0.1; % vert. space between subplots
lefs = 0.15; % very left of figure
rigs = 0.12; % very right of figure
tops = 0.02; % top of figure
bots = 0.16; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

pcolor(xVecCTD, zVec, O2(:,theIndexCTD)); shading interp
hold on
contour(xVecCTD, zVec, sig0(:,theIndexCTD), v1, 'color', [.5 .5 .5])
map = brewermap(14, 'RdBu');
colormap(map)
caxis([100, 250])
ylabel('Depth (m)', 'FontSize', FS1, 'fontWeight', 'bold')
set(gca, 'ygrid', 'on');
set(gca, 'ydir', 'reverse');
set(gca, 'xtick', XTICKS)
set(gca, 'tickdir', 'out');
xlim(XLIM)
ylim([0 zMax])
cb = colorbar;
text(XLIM(1)+2,zMax-.1*zMax, 'O_2 (\mumol Kg^{-1})  (I)', 'horizontalAlignment', 'left', 'fontSize', FS1, 'fontWeight', 'bold', 'color',[1 1 1]*0)

adjust_space
pause(1)

cpos = get(cb, 'pos');
cpos(1) = cpos(1) - cbar_offset;
cpos(2) = cpos(2)+(cpos(4)-cpos(4)*vredux)/2;
cpos(4) = cpos(4)*vredux;
cpos(3) = cpos(3)*hredux;
set(cb, 'pos', cpos)
set(gcf, 'renderer', 'painters')
print(gcf, '-dpng', '-r300', 'M327_O2_noplot.png')  


pause(1)

!montage M327_T_noplot.png M327_S_noplot.png M327_O2_noplot.png M327_CHL_noplot.png M327_TRY_noplot.png M327_PHE_noplot.png M327_CDOM_noplot.png M327_BB_noplot.png -tile 2x4  -geometry 740x390+1+1  M327_preswot.png 