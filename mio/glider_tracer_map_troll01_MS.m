%% THIS SCRIPT GENERATES THE FIGURE FOUND IN MANUSCRIPT ON PAH
% -> Run glider_proj_troll01 first!

%% some plot properties
cbar_offset = -0.1;
FONTSIZE = 10;
vredux = .8;
hredux = .7;


%% By-pass existing calibration   <------------ ***
% (ideally, netcdf files should be re-generated with new calib.)
NAP_calib = 0.3154; % N=5
NAP_blank = 0.0684;
PHE_calib = 3.0869;
PHE_blank = -0.0177;
%I = find(TRYru>.14); % <---- Try to flag bad counts See how I can do it better!
%TRYru(I) = NaN;
PHE = ( PHEru - PHE_blank)./PHE_calib;
NAP = ( TRYru - NAP_blank)./NAP_calib;

%% Manual removal of problematic counts in NAP
NAP(:,36) = NaN;
NAP(75:300,92) = NaN;
NAP(1:40,72) = NaN;
NAP(1:100,74) = NaN;


%% Deal with map limits
region = 'rockall'; % for the bathym
LIMS = [3.18 4.4 60.6 61]; % lat/lon
lonLims = LIMS(1:2);
latLims = LIMS(3:4);

%% bathymetry
path = ['~/data/matlab_bathym/' region '.mat']; % 30-sec.
load(path)
I=find(lat<latLims(2) & lat>latLims(1));
J=find(lon<lonLims(2) & lon>lonLims(1));
latitude=lat(I);
longitude=lon(J);
bathy=z(I,J);       

%% plot map
close all
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 2; % no. subplot row
dx = 0.0 ; % horiz. space between subplots
dy = 0.0; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.15; % very right of figure
tops = 0.01; % top of figure
bots = 0.08; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %
FS = 10;
%V = 100:1:300;

figure(1);
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 17 20])

% -> Naphs
tracer = nanmean(NAP(1:10, :));
subplot(211)
m_proj('mercator','long',lonLims,'lat',latLims);
hold on
load('PuBuGn_modif.mat');
%colormap(PuBuGn_modif);
colormap(brewermap(8,'PuBuGn'))
%colormap(parula);
V=[0:100:1000];
[cc, HH] = m_contour(longitude,latitude,-bathy, [0:500:3000], 'color', 'k');
m_gshhs_h('patch',[1 .9333 .6667]); %coastlines (Beige)
clabel(cc, HH)   
ylabel('Latitude', 'FontSize', FS, 'fontweight', 'bold')

[x,y] = m_ll2xy(lonVec, latVec);
h=scatter(x,y,25,tracer,'filled');
caxis([0.05, .16])
c = colorbar;
ylabel(c,'Naph-like (ug L^{-1})', 'FontSize', FS, 'fontweight', 'bold')
m_grid('box','fancy')


% Troll ABC
trollA = [60.6563, 3.7270];
trollB = [60.7752, 3.4922];
trollC = [60.8883, 3.6090];
SEn = [60.7525, 3.4373];
CP = [60.7673, 3.6082];
SEq = [60.8302, 3.5699];
HT = [60.7109, 3.4040];
m_line(trollA(2), trollA(1), 'marker', 'p', 'color', 'm', 'markerFaceColor', ...
       'm', 'markersize',14);
m_line(trollB(2), trollB(1), 'marker', 'p', 'color', 'm', 'markerFaceColor', ...
       'm', 'markersize',14);
m_line(trollC(2), trollC(1), 'marker', 'p', 'color', 'm', 'markerFaceColor', ...
       'm', 'markersize',14);
m_line(SEn(2), SEn(1), 'marker', 'p', 'color', 'b', 'markerFaceColor', ...
       'b', 'markersize',14);
m_line(CP(2), CP(1), 'marker', 'p', 'color', 'b', 'markerFaceColor', ...
       'b', 'markersize',14);
m_line(SEq(2), SEq(1), 'marker', 'p', 'color', 'b', 'markerFaceColor', ...
       'b', 'markersize',14);
m_line(HT(2), HT(1), 'marker', 'p', 'color', 'c', 'markerFaceColor', ...
       'c', 'markersize',14);
m_text(trollA(2)+.01, trollA(1), 'A', 'color', 'm', 'horizontalAlignment', 'left', ...
       'verticalAlignment', 'top')
m_text(trollB(2)+.01, trollB(1), 'B', 'color', 'm', 'horizontalAlignment', 'left', ...
       'verticalAlignment', 'top')
m_text(trollC(2)+.01, trollC(1), 'C', 'color', 'm', 'horizontalAlignment', 'left', ...
       'verticalAlignment', 'top')
m_text(SEn(2)+.01, SEn(1), 'SEn', 'color', 'b', 'horizontalAlignment', 'left', ...
       'verticalAlignment', 'top')
m_text(CP(2)+.01, CP(1), 'CP', 'color', 'b', 'horizontalAlignment', 'left', ...
       'verticalAlignment', 'top')
m_text(SEq(2)+.01, SEq(1), 'SEq', 'color', 'b', 'horizontalAlignment', 'left', ...
       'verticalAlignment', 'top')
m_text(HT(2)+.01, HT(1), 'HT', 'color', 'c', 'horizontalAlignment', 'left', ...
       'verticalAlignment', 'top', 'fontWeight', 'bold')
adjust_space

pause(1)
%keyboard
cpos = get(c, 'pos');
cpos(1) = cpos(1) - cbar_offset;
cpos(2) = cpos(2)+(cpos(4)-cpos(4)*vredux)/2;
cpos(4) = cpos(4)*vredux;
cpos(3) = cpos(3)*hredux;
set(c, 'pos', cpos)



% -> Phes
tracer = nanmean(PHE(1:10, :))*1000;
subplot(212)
m_proj('mercator','long',lonLims,'lat',latLims);
hold on
V=[0:100:1000];
[cc, HH] = m_contour(longitude,latitude,-bathy, [0:500:3000], 'color', 'k');
m_gshhs_h('patch',[1 .9333 .6667]); %coastlines (Beige)
clabel(cc, HH)   
ylabel('Latitude', 'FontSize', FS, 'fontweight', 'bold')
xlabel('Longitude', 'FontSize', FS, 'fontweight', 'bold')

[x,y] = m_ll2xy(lonVec, latVec);
h=scatter(x,y,25,tracer,'filled');
caxis([22, 28])
c = colorbar;
ylabel(c,'Phe-like (ng L^{-1})', 'FontSize', FS, 'fontweight', 'bold')
m_grid('box','fancy')


% Troll ABC
trollA = [60.6563, 3.7270];
trollB = [60.7752, 3.4922];
trollC = [60.8883, 3.6090];
SEn = [60.7525, 3.4373];
CP = [60.7673, 3.6082];
SEq = [60.8302, 3.5699];
HT = [60.7109, 3.4040];
m_line(trollA(2), trollA(1), 'marker', 'p', 'color', 'm', 'markerFaceColor', ...
       'm', 'markersize',14);
m_line(trollB(2), trollB(1), 'marker', 'p', 'color', 'm', 'markerFaceColor', ...
       'm', 'markersize',14);
m_line(trollC(2), trollC(1), 'marker', 'p', 'color', 'm', 'markerFaceColor', ...
       'm', 'markersize',14);
m_line(SEn(2), SEn(1), 'marker', 'p', 'color', 'b', 'markerFaceColor', ...
       'b', 'markersize',14);
m_line(CP(2), CP(1), 'marker', 'p', 'color', 'b', 'markerFaceColor', ...
       'b', 'markersize',14);
m_line(SEq(2), SEq(1), 'marker', 'p', 'color', 'b', 'markerFaceColor', ...
       'b', 'markersize',14);
m_line(HT(2), HT(1), 'marker', 'p', 'color', 'c', 'markerFaceColor', ...
       'c', 'markersize',14);
adjust_space
pause(1)

cpos = get(c, 'pos');
cpos(1) = cpos(1) - cbar_offset;
cpos(2) = cpos(2)+(cpos(4)-cpos(4)*vredux)/2;
cpos(4) = cpos(4)*vredux;
cpos(3) = cpos(3)*hredux;
set(c, 'pos', cpos)


set(gcf, 'renderer', 'painters')
print('-dpng', '-r300',  'surface_PAH_troll.png')

