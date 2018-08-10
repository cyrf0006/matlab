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
lefs = 0.05; % very left of figure
rigs = 0.05; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

%% plotting info
%figure dimension
paperwidth = 25;%cm
paperheight = 20;%cm
FS = 15;

% param for colorbar
cbar_width = 0.02;
cbar_offset = 0.01; % colorbar offset from figure
offset2 = 0.1; 


%% Get global bathym    
path = ['~/data/matlab_bathym/rockall.mat']; % 30-sec.
lims = [40 70 -35 10];
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
%decimFactor = 4;
decimFactor = 40;
lat1m = lat(1:decimFactor:end);
lon1m = lon(1:decimFactor:end);
z1m = z(1:decimFactor:end, 1:decimFactor:end);

%% Zoom map
lon_min=-20;
lon_max=-11;
lat_min=54;
lat_max=59;

I=find(lat<lat_max & lat>lat_min);
J=find(lon<lon_max & lon>lon_min);
latitude=lat(I);
longitude=lon(J);
bathy=z(I,J);

%% AVISO tides
tide = load('~/research/CimaStuff/tide/tide_Oct_2012.mat');
lonVecTide = tide.lon_range;
latVecTide = tide.lat_range;
timeVecTide = tide.time;
%etaTide = tide.slev;
origin = datenum(1950, 1, 1);
timeVecTide = origin + timeVecTide;

[Y, I] = min(abs(timeVecTide-datenum(2012,10,16,18,0,0)));
[Y, I] = min(abs(timeVecTide-datenum(2012,10,17,03,0,0)));
vMat = squeeze(tide.vv(:,:,I));
uMat = squeeze(tide.uu(:,:,I));

[uMat2, vMat2] = rotate_vecd(uMat, vMat, 62);
vMat = uMat2;


%% ------- Plot main figure ------ %%
%h = figure('visible', 'off');
h = figure(1);
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 paperwidth paperheight])
m_proj('mercator','long',[min(longitude) max(longitude)],'lat',[min(latitude) max(latitude)]);
%load ~/PhD/bathym/gebco64
%colormap(gebco);
load BR_symetric
colormap(BR_symetric);

V=[0:100:4000];
%m_pcolor(longitude,latitude,abs(bathy)); shading flat;
m_contourf(lonVecTide,latVecTide,vMat, 50, 'lineStyle', 'none');
%shading flat
hold on
[HH, HH] = m_contour(longitude,latitude,abs(bathy), [0:500:2500], 'color', 'k');
m_gshhs_h('patch',[1 .9333 .6667]); %coastlines (Beige)                               
xlabel('Longitude', 'FontSize', FS, 'fontweight', 'bold')
%ylabel('Latitude', 'FontSize', FS, 'fontweight', 'bold')
%caxis([min(V) max(V)])
set(gca, 'fontsize', FS)
caxis([-25 25])

% Rectangle around study region (zoom++ )
lon_min=-15.85;
lon_max=-15.75;
lat_min=55.45;
lat_max=55.5333;

[lomi, lami] = m_ll2xy(lon_min,lat_min);
[loma, lama] = m_ll2xy(lon_max,lat_max);
rectangle('Position', [lomi lami loma-lomi lama-lami],'linewidth', 1, 'edgecolor', 'r') 
m_grid('box','fancy', 'yticklabel', [])

adjust_space

%% get bathym transect
m_line([-16.1120, -15.3720],[55.820, 55.023],'color','m','linewi',2,'linestyle', '--'); 
[range,ln,lt]=m_lldist([-16.1120, -15.3720],[55.820, 55.023], 5000); %'m'



set(gcf, 'renderer', 'painters'); % vectorial figure
print('-depsc', 'RockallMapAviso_plain.eps')


keyboard

