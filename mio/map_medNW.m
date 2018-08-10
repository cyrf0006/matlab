clear all
close all
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.08; % very left of figure
rigs = 0.05; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

% plotting info
%figure dimension
paperwidth = 24;%cm
paperheight = 15;%cm
FS = 14;

% param for colorbar
cbar_width = 0.02;
cbar_offset = 0.01; % colorbar offset from figure
offset2 = 0.1; 


% ZOOM LIMS
lonLims = [0 12];
latLims = [40 45];   

% bathymetry
path = ['~/data/matlab_bathym/MedWest.mat']; % 30-sec.
load(path)
I=find(lat<latLims(2) & lat>latLims(1));
J=find(lon<lonLims(2) & lon>lonLims(1));
latitude=lat(I);
longitude=lon(J);
bathy=z(I,J);


%% Main pannel %%
h = figure('visible', 'off');
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 paperwidth paperheight])
m_proj('mercator', 'long',lonLims,'lat',latLims);
%m_grid('box','fancy')
hold on
V=[0:100:4000];
%[HH, HH] = m_contour(longitude,latitude,-bathy, [0:500:4000], 'color', 'k');
m_gshhs_f('patch',[1 1 1]*.9);
xlabel('Longitude', 'FontSize', FS, 'fontweight', 'bold')
ylabel('Latitude', 'FontSize', FS, 'fontweight', 'bold')
m_grid('box','fancy')

% Rectangle
[lomi, lami] = m_ll2xy(6,42.333);
[loma, lama] = m_ll2xy(9,44);
rectangle('Position', [lomi lami loma-lomi lama-lami],'linewidth', 1, 'edgecolor', 'r') 

% add political boundaries
M=m_shaperead('/home/cyrf0006/data/naturalEarth/ne_10m_admin_0_boundary_lines_land'); 
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

adjust_space




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

%save figure
print('-dpng', '-r300', 'map_medNW.png')
%set(gcf, 'renderer', 'painters'); % vectorial figure
%print('-depsc', 'RockallMap.eps')


keyboard

