clear all
close all

%% plotting info
%figure dimension
paperwidth = 12;%cm
paperheight = 15;%cm
FS = 14;
load gebco64

%% mission parameters
wp_recov = [39, 3];
track713_0 = [36.039826, 3.923135];
track713_1 = [39.959253, 2.683943];
orig_1 = [40.002187929906384, 2.701594620126058];
orig_2 = [39.91631506674032, 2.666313445481265];

end_1 = [36.08278979266785, 3.939869792742548];
end_2 = [35.99685956567776, 3.9064183788794073];


SW1 = [39.0, 2.97];
SWf = [37.0, 3.60];
NE1 = [39.0, 3.03];
NEf = [37.0, 3.65];

% How I get this with python tool:
% $$$ In [1]: import geo_tools as g
% $$$ In [2]: lat1 = 39.959253
% $$$ In [3]: lon1 = 2.683943
% $$$ In [4]: distance=5
% $$$ In [5]: bearing = g.azimuth(2.683943, 39.959253, 3.923135, 36.039826)
% $$$ In [8]: g.get_azimuth_point(lat1, lon1, bearing-90, distance)
% $$$ Out[8]: Point(40.002187929906384, 2.701594620126058, 0.0)
% $$$ In [9]: g.get_azimuth_point(lat1, lon1, bearing+90, distance)
% $$$ Out[9]: Point(39.91631506674032, 2.666313445481265, 0.0)


% param for colorbar
cbar_width = 0.02;
cbar_offset = 0.01; % colorbar offset from figure
offset2 = 0.1; 

lonmin=2.5;
lonmax=4;
latmin=36.5;
latmax=39.5;
% ZOOM LIMS
lonLims = [lonmin lonmax];
latLims = [latmin latmax];   

% bathymetry
path = ['~/Data/matlab_bathym/westMed.mat']; % 30-sec.
lims = [35 45 -6 9];
fid = fopen(path);
if fid == -1 % doesnt exist yet!
    disp('30-sec. region doesn''t exist yet!')
    disp(' Extracting from GEBCO, may take some time...')
    [lat lon z] = getGebco('~/Data/GEBCO/GEBCO_08.nc', 'z', lims);
    disp('   -> done!')
    save(path, 'lat', 'lon', 'z');
else
    load(path)
end
load(path)
I=find(lat<latLims(2) & lat>latLims(1));
J=find(lon<lonLims(2) & lon>lonLims(1));
latitude=lat(I);
longitude=lon(J);
bathy=z(I,J);


%% Main pannel %%
h = figure('visible', 'on');
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 paperwidth paperheight])
m_proj('mercator', 'long',lonLims,'lat',latLims);
%m_grid('box','fancy')
hold on
V=[0:100:4000];
[HH, HH] = m_contourf(longitude,latitude,-bathy, [0:100:4000], 'lineStyle', 'none');
[HH, HH] = m_contour(longitude,latitude,-bathy, [0:500:4000], 'color', 'k');
m_gshhs_h('patch',[1 1 1]*.9);
xlabel('Longitude', 'FontSize', FS, 'fontweight', 'bold')
ylabel('Latitude', 'FontSize', FS, 'fontweight', 'bold')
m_grid('box','fancy')
colormap(gebco)


% add political boundaries
M=m_shaperead('/home/cyrf0006/Data/naturalEarth/ne_10m_admin_0_boundary_lines_land'); 
for k=1:length(M.ncst), 
  m_line(M.ncst{k}(:,1),M.ncst{k}(:,2), 'color', [1 1 1]*.5); 
end; 

% add cities
m_plot(7.26, 43.7, 'p', 'color', 'k', 'markerFaceColor', 'k', 'markerSize', 8);
m_text(7.24, 43.71, 'Nice', 'color', 'k', 'fontSize', 10, 'fontWeight', ...
       'bold', 'horizontalAlignment', 'right', 'verticalAlignment', 'bottom');
m_plot(8.76, 42.55, 'p', 'color', 'k', 'markerFaceColor', 'k', 'markerSize', 8);
m_text(8.76, 42.55, 'Calvi', 'color', 'k', 'fontSize', 10, 'fontWeight', ...
       'bold', 'horizontalAlignment', 'left', 'verticalAlignment', 'top');
m_text(6.1, 42.4, 'Study Region', 'color', 'r', 'fontSize', 10, 'fontWeight', ...
       'bold', 'horizontalAlignment', 'left', 'verticalAlignment', 'bottom');
m_text(1, 41.75, 'Spain', 'color', 'k', 'fontSize', 14, 'horizontalAlignment', 'left', 'verticalAlignment', 'middle');
m_text(4.5, 44.5, 'France', 'color', 'k', 'fontSize', 14, 'horizontalAlignment', 'left', 'verticalAlignment', 'middle');
m_text(10.1, 44.5, 'Italy', 'color', 'k', 'fontSize', 14, 'horizontalAlignment', 'left', 'verticalAlignment', 'middle');

% Plot satellite tracks
plot_Sentinel('./S3A_C0030_T0244_LON-002+006_LAT+36+40.txt.compatible',lonmin,lonmax,latmin,latmax,'r--')
plot_Sentinel('./S3A_C0030_T0713_LON-002+006_LAT+36+40.txt.compatible',lonmin,lonmax,latmin,latmax,'b--')

m_plot(wp_recov(2), wp_recov(1), 'gp', 'markerSize', 12, 'markerFaceColor', 'g')

m_plot([SW1(2), SWf(2)], [SW1(1), SWf(1)], '--k', 'lineWidth', 2)
m_plot([NE1(2), NEf(2)], [NE1(1), NEf(1)], '--k', 'lineWidth', 2)
[R, LN, LT] = m_lldist([SW1(2), SWf(2)], [SW1(1), SWf(1)], 5);
m_plot(LN, LT, 'om', 'markerFaceColor', 'm')
for i=1:length(LN)
    m_text(LN(i), LT(i), sprintf('SW%d  ', i), 'fontWeight', ...
           'bold', 'fontSize', 10, 'horizontalAlignment', 'right')
end
disp('SouthWest track [Lat, Lon]:')
[LT', LN']

[R, LN, LT] = m_lldist([NE1(2), NEf(2)], [NE1(1), NEf(1)], 5);
m_plot(LN, LT, 'or', 'markerFaceColor', 'r')
for i=1:length(LN)
    m_text(LN(i), LT(i), sprintf('  NE%d', i), 'fontWeight', ...
           'bold', 'fontSize', 10, 'horizontalAlignment', 'left') 
end
disp('NorthEast track [Lat, Lon]:')
[LT', LN']
% Draw scale km
m_plot(LN(1:2)+.4, LT(1:2), '-k', 'linewidth', 3)
m_text(mean(LN(1:2))+.5, mean(LT(1:2)), '45.8km', 'fontWeight', 'bold','fontSize', 10, 'horizontalAlignment', 'left')


disp('Distance Between tracks (km):')
m_lldist([NE1(2), SW1(2)], [NE1(1), SW1(1)])  

disp('Coordinates of the recovery point [Lat, Lon]')
wp_recov


%% Inset %%
lon_min=-5;
lon_max=20;
lat_min=35;
lat_max=48;

a2 = axes('position',[0.08 0.65 0.26 0.25]) ; % inset
m_proj('mercator', 'long',[lon_min lon_max],'lat',[lat_min lat_max]);
hold on
m_gshhs_i('patch',[1 1 1]*.9);
caxis([min(V) max(V)])
m_grid('box','fancy', 'yticklabels',[], 'xticklabels',[])

% Rectangle
[lomi, lami] = m_ll2xy(lonLims(1),latLims(1));
[loma, lama] = m_ll2xy(lonLims(2),latLims(2));

rectangle('Position', [lomi lami loma-lomi lama-lami],'linewidth', 1, 'edgecolor', 'r') 
set(gca, 'fontSize', FS)



print('-dpng', '-r300', 'map_bioswot.png')








