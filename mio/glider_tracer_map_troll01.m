function glider_tracer_map_troll01(region, latVec, lonVec, LIMS, tracer, CLIM, TITLE, outfile)

% to be run after glider_contours2.m
% glider_tracer_map_troll01('rockall', latVec, lonVec, [3.18 4.4 60.6 61],nanmean(PHE(1:10, :))*1000, [20 35], 'Phe (ng/L)', 'Phe_troll_surface.png')
% glider_tracer_map_troll01('rockall', latVec, lonVec, [3.18 4.4 60.6 61],nanmean(PHE_ru(1:10, :)), [.04 .08], 'Phe (RU)', 'Phe_troll_surface_RU.png')
% glider_tracer_map_troll01('rockall', latVec, lonVec, [3.18 4.4 60.6 61],nanmean(TRY_ru(1:10, :)), [.06 .14], 'Naph (RU)', 'Naph_troll_surface_RU.png')

%% Deal with varargin (map limits)
lonLims = LIMS(1:2);
latLims = LIMS(3:4);
% $$$ lonLims = [min(lonVec)-.5 max(lonVec)+1];
% $$$ latLims = [min(latVec)-.7 max(latVec)+1];        

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
nrow = 1; % no. subplot row
dx = 0.0 ; % horiz. space between subplots
dy = 0.0; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.13; % very right of figure
tops = 0.05; % top of figure
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
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 18 12])
m_proj('mercator','long',lonLims,'lat',latLims);
hold on
colormap(jet);
V=[0:100:1000];
%m_pcolor(longitude,latitude,-bathy); shading flat;
[cc, HH] = m_contour(longitude,latitude,-bathy, [0:500:3000], 'color', 'k');
m_gshhs_h('patch',[1 .9333 .6667]); %coastlines (Beige)
clabel(cc, HH)   
xlabel('Longitude', 'FontSize', FS, 'fontweight', 'bold')
ylabel('Latitude', 'FontSize', FS, 'fontweight', 'bold')


[x,y] = m_ll2xy(lonVec, latVec);
h=scatter(x,y,25,tracer,'filled');
caxis(CLIM)
c = colorbar;

%m_plot(oilPatchLon, oilPatchLat, '--', 'color', 'g', 'lineWidth', 2);

title(TITLE)
m_grid('box','fancy')



%% add platforms
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

print('-dpng', '-r300',  outfile)

