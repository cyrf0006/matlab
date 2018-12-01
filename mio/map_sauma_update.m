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
lefs = 0.06; % very left of figure
rigs = 0.03; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

%% plotting info
%figure dimension
paperwidth = 18;%cm
paperheight = 20;%cm
FS = 12;
FS = 12;
FS2 = 10;
% param for colorbar
cbar_width = 0.02;
cbar_offset = 0.01; % colorbar offset from figure
offset2 = 0.1; 

% Troll ABC
trollA = [60.6563, 3.7270];
trollB = [60.7752, 3.4922];
trollC = [60.8883, 3.6090];
SEn = [60.7525, 3.4373];
CP = [60.7673, 3.6082];
SEq = [60.8302, 3.5699];

WP = ...
    [ 60+45.0/60 3+12.0/60
      60+55.0/60 4+22.0/60
      60+46.8/60 3+24.6/60
      60+41.0/60 3+47.0/60
      60+52.0/60 3+35.2/60];
WP_names =  {'L1', 'L2', 'T1', 'T2' 'T3'};

lonVec = load('./gliderLon.txt');
latVec = load('./gliderLat.txt');

%% Get global bathym    
path = ['~/data/matlab_bathym/troll.mat']; % 30-sec.
path = ['~/data/matlab_bathym/rockall.mat']; % 30-sec.
%lims = [50 65  -10 15];
lims = [50 80 -35 20];

fid = fopen(path);
if fid == -1 % doesnt exist yet!
    disp('30-sec. region doesn''t exist yet!')
    disp(' Extracting from GEBCO, may take some time...')
    [lat lon z] = getGebco('~/data/GEBCO/GEBCO_08.nc', 'z', lims);
    disp('   -> done!')
    save(path, 'lat', 'lon', 'z');
else
    load(path)
end

% Lower res.
decimFactor = 40;
lat1m = lat(1:decimFactor:end);
lon1m = lon(1:decimFactor:end);
z1m = z(1:decimFactor:end, 1:decimFactor:end);



%% Zoom map
lon_min=-6;
lon_max=15;
lat_min=54;
lat_max=63;

I=find(lat<lat_max & lat>lat_min);
J=find(lon<lon_max & lon>lon_min);
latitude=lat(I);
longitude=lon(J);
bathy=z(I,J);
I = find(bathy>0);
bathy(I) = 0;

h = figure('visible', 'off');
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 paperwidth paperheight])
m_proj('mercator','long',[min(longitude) max(longitude)],'lat',[min(latitude) max(latitude)]);
%m_grid('box','fancy')
hold on
load gebco64
colormap(gebco);
V=[0:50:1000];
m_pcolor(longitude,latitude,abs(bathy)); shading flat;
[HH, HH] = m_contour(longitude,latitude,abs(bathy), [0:100:1000], 'color', 'k');
m_gshhs_h('patch',[1 .9333 .6667]); %coastlines (Beige)                               
xlabel('Longitude', 'FontSize', FS, 'fontweight', 'bold')
ylabel('Latitude', 'FontSize', FS, 'fontweight', 'bold')
caxis([min(V) max(V)])

% Rectangle around study region (zoom++ )
lon_min=3.15;
lon_max=4.4;
lat_min=60.5;
lat_max=61;

[lomi, lami] = m_ll2xy(lon_min,lat_min);
[loma, lama] = m_ll2xy(lon_max,lat_max);
rectangle('Position', [lomi lami loma-lomi lama-lami],'linewidth', 1, 'edgecolor', 'r') 

% adjust_these lines;  
m_line([lon_min, -5.6],[lat_max, 58.35],'color','r','linewi',1,'linestyle', '-');        
m_line([lon_max, 3.5],[lat_min, 54.25],'color','r','linewi',1,'linestyle', '-');  
m_grid('box','fancy')


m_line(trollA(2), trollA(1), 'marker', 'p', 'color', 'b', 'markerFaceColor', ...
       'b', 'markersize',14);
m_line(trollB(2), trollB(1), 'marker', 'p', 'color', 'b', 'markerFaceColor', ...
       'b', 'markersize',14);
m_line(trollC(2), trollC(1), 'marker', 'p', 'color', 'b', 'markerFaceColor', ...
       'b', 'markersize',14);
m_line(SEn(2), SEn(1), 'marker', 'p', 'color', 'b', 'markerFaceColor', ...
       'b', 'markersize',14);
m_line(CP(2), CP(1), 'marker', 'p', 'color', 'b', 'markerFaceColor', ...
       'b', 'markersize',14);
m_line(SEq(2), SEq(1), 'marker', 'p', 'color', 'b', 'markerFaceColor', ...
       'b', 'markersize',14);

adjust_space


%% Large Area map (inset)
lon_min=-15;
lon_max=20;
lat_min=50;
lat_max=70;

I=find(lat1m<lat_max & lat1m>lat_min);
J=find(lon1m<lon_max & lon1m>lon_min);
latitude=lat1m(I);
longitude=lon1m(J);
bathy=z1m(I,J);
I = find(bathy>0);
bathy(I) = 0;

a2 = axes('position',[0.1 0.6 0.33 0.33]) ; % inset
m_proj('mercator', 'long',[lon_min lon_max],'lat',[lat_min lat_max]);
hold on
[HH, HH] = m_contour(longitude,latitude,abs(bathy), [0:1000:4000], 'color', 'k');
m_gshhs_i('patch',[1 .9333 .6667]); %coastlines (Beige)                               
caxis([min(V) max(V)])
m_grid('box','fancy', 'yticklabels',[], 'xticklabels',[])

% Rectangle
lon_min=-6;
lon_max=15;
lat_min=54;
lat_max=63;
[lomi, lami] = m_ll2xy(lon_min,lat_min);
[loma, lama] = m_ll2xy(lon_max,lat_max);

rectangle('Position', [lomi lami loma-lomi lama-lami],'linewidth', 1, 'edgecolor', 'r') 


% --- Colorbar --- %
cbar_offset = 0;%0.04; 
offset2 =0;% 0.02; % offset between heigth of colorbar 
c = colorbar('location','southoutside');
cb_pos = get(c, 'position');
cb_back = cb_pos;
cb_pos(4) = .02;
cb_pos(1) = .2;%cb_pos(1)+.1*cb_pos(3);
cb_pos(3) = .32;
cb_pos(2) = .52;%cb_pos(2)+cbar_offset;
set(c, 'pos', cb_pos);
set(c, 'fontsize', FS2, 'fontweight', 'bold')

ti = ylabel(c,'Depth (m)', 'FontSize', FS, 'fontweight', 'bold');
ti_pos = get(ti, 'position');
set(ti, 'Rotation',0.0);     
clim = get(gca, 'clim');
set(ti, 'position', [(clim(2)-clim(1))./2  2 ti_pos(3)]); 
set(gca, 'fontSize', FS)



%% Study region Map (zoom ++)
lon_min=3.15;
lon_max=4.4;
lat_min=60.5;
lat_max=61;

I=find(lat<lat_max & lat>lat_min);
J=find(lon<lon_max & lon>lon_min);
latitude=lat(I);
longitude=lon(J);
bathy=z(I,J);
I = find(bathy>0);
bathy(I) = 0;
V = 0:100:400;
%a2 = axes('position',[0.55 0.12 0.4 0.3]) ; % inset
a2 = axes('position',[0.1 0.1 0.5 0.4]) ; % inset
m_proj('mercator', 'long',[lon_min lon_max],'lat',[lat_min lat_max]);
[HH, HH] = m_contour(longitude,latitude,abs(bathy), V, 'color', 'k');
hold on
m_gshhs_f('patch',[1 .9333 .6667]); %coastlines (Beige)                               
caxis([min(V) max(V)])
set(HH, 'color', 'k', 'ShowText','on')%, 'LabelSpacing', 800)

%% add track
for i = 1:length(lonVec)
    m_plot(lonVec(i),latVec(i), '.', 'color', 'r', 'markerSize', 10);
end
    m_plot(lonVec(end),latVec(end), '.', 'color', 'b', 'markerSize', 14);

%% add platforms
m_line(trollA(2), trollA(1), 'marker', 'p', 'color', 'b', 'markerFaceColor', ...
       'b', 'markersize',14);
m_line(trollB(2), trollB(1), 'marker', 'p', 'color', 'b', 'markerFaceColor', ...
       'b', 'markersize',14);
m_line(trollC(2), trollC(1), 'marker', 'p', 'color', 'b', 'markerFaceColor', ...
       'b', 'markersize',14);
m_line(SEn(2), SEn(1), 'marker', 'p', 'color', 'b', 'markerFaceColor', ...
       'b', 'markersize',14);
m_line(CP(2), CP(1), 'marker', 'p', 'color', 'b', 'markerFaceColor', ...
       'b', 'markersize',14);
m_line(SEq(2), SEq(1), 'marker', 'p', 'color', 'b', 'markerFaceColor', ...
       'b', 'markersize',14);
% $$$ m_text(trollA(2)+.01, trollA(1), 'A', 'color', 'm', 'horizontalAlignment', 'left', ...
% $$$        'verticalAlignment', 'top')
% $$$ m_text(trollB(2)+.01, trollB(1), 'B', 'color', 'm', 'horizontalAlignment', 'left', ...
% $$$        'verticalAlignment', 'top')
% $$$ m_text(trollC(2)+.01, trollC(1), 'C', 'color', 'm', 'horizontalAlignment', 'left', ...
% $$$        'verticalAlignment', 'top')
% $$$ % $$$ m_text(SEn(2)+.01, SEn(1), 'SEn', 'color', 'b', 'horizontalAlignment', 'left', ...
% $$$        'verticalAlignment', 'top')
% $$$ m_text(CP(2)+.01, CP(1), 'CP', 'color', 'b', 'horizontalAlignment', 'left', ...
% $$$        'verticalAlignment', 'top')
% $$$ m_text(SEq(2)+.01, SEq(1), 'SEq', 'color', 'b', 'horizontalAlignment', 'left', ...
% $$$        'verticalAlignment', 'top')

%% add WPs
m_plot([WP(end-4:end,2); WP(end-2,2)],[WP(end-4:end,1); WP(end-2,1)], '--k', 'lineWidth', 2)
for i = 1:size(WP,1)
    m_plot(WP(i,2),WP(i,1), '.', 'color', [.3 .3 .3], 'markerSize', 14);
end

m_text(WP(1,2),WP(1,1), WP_names(1), 'color', [.3 .3 .3], 'verticalAlignment', 'bottom', 'horizontalAlignment', 'center', 'fontSize', 14, 'fontweight', 'bold');
m_text(WP(2,2),WP(2,1), WP_names(2), 'color', [.3 .3 .3], 'verticalAlignment', 'top', 'horizontal', 'right', 'fontSize', 14, 'fontweight', 'bold');
m_text(WP(3,2),WP(3,1), WP_names(3), 'color', [.3 .3 .3], 'verticalAlignment', 'bottom', 'horizontal', 'right', 'fontSize', 14, 'fontweight', 'bold');
m_text(WP(4,2),WP(4,1), WP_names(4), 'color', [.3 .3 .3], 'verticalAlignment', 'middle', 'horizontalAlignment', 'left', 'fontSize', 14, 'fontweight', 'bold');
m_text(WP(5,2),WP(5,1), WP_names(5), 'color', [.3 .3 .3], 'verticalAlignment', 'bottom', 'horizontal', 'right', 'fontSize', 14, 'fontweight', 'bold');

m_grid('box','fancy', 'yticklabels',[], 'xticklabels',[])



%save figure
print('-dpng', '-r300', 'mapTroll_update.png')
%set(gcf, 'renderer', 'painters'); % vectorial figure
%print('-depsc', 'mapTroll.eps')

