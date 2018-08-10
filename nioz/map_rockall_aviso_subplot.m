clear all
close all
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 2; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.02 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.05; % very left of figure
rigs = 0.05; % very right of figure
tops = 0.1; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

%% plotting info
%figure dimension
paperwidth = 40;%cm
paperheight = 25;%cm
FS = 15;

% param for colorbar
cbar_width = 0.02;
cbar_offset = 0.01; % colorbar offset from figure
offset2 = 0.1; 

% time of the tide
t0 = datenum(2012,10,17,03,0,0);

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

%% Get local bathym (multibeam)
path2 = ['~/research/NIOZ/RockallBank/bathymetry/XYZHERMESALL_ROCKALL.CSV'];
lims = [40 70 -35 10];
CWC = load(path2);

latCWC = CWC(:,1);
lonCWC = CWC(:,2);
zCWC = CWC(:,3);
%reduce size (zoom ++)
lon_min=-15.85;
lon_max=-15.6;
lat_min=55.45;
lat_max=55.55;
I=find(latCWC<=lat_max & latCWC>=lat_min & lonCWC<=lon_max & lonCWC>=lon_min);
latCWC = latCWC(I);
lonCWC = lonCWC(I);
zCWC = zCWC(I);


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

[Y, I] = min(abs(timeVecTide-(t0-.5)));
time1 = timeVecTide(I);
vMat = squeeze(tide.vv(:,:,I));
uMat = squeeze(tide.uu(:,:,I));

[uMat2, vMat2] = rotate_vecd(uMat, vMat, 62);
vMat = uMat2;


%% ------- Plot main figure ------ %%
%h = figure('visible', 'off');
h = figure(1);
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 paperwidth paperheight])

%% SUBPLOT 1
s = subplot(121);
m_proj('mercator','long',[min(longitude) max(longitude)],'lat',[min(latitude) max(latitude)]);
load BR_symetric
colormap(BR_symetric);

V=[0:100:4000];
m_contourf(lonVecTide,latVecTide,vMat/100, 50, 'lineStyle', 'none');
hold on
[HH, HH] = m_contour(longitude,latitude,abs(bathy), [0:500:2500], 'color', 'k');
m_gshhs_h('patch',[1 .9333 .6667]); %coastlines (Beige)                               
xlabel('Longitude', 'FontSize', FS, 'fontweight', 'bold')
ylabel('Latitude', 'FontSize', FS, 'fontweight', 'bold')
set(gca, 'fontsize', FS)
caxis([-.25 .25])

% Rectangle around study region (zoom++ )
lon_min=-15.85;
lon_max=-15.6;
lat_min=55.45;
lat_max=55.55;

[lomi, lami] = m_ll2xy(lon_min,lat_min);
[loma, lama] = m_ll2xy(lon_max,lat_max);
rectangle('Position', [lomi lami loma-lomi lama-lami],'linewidth', 1, 'edgecolor', 'r') 

% few annotations
m_text(-14.75, 56.6667, 'Rockall Bank', 'rotation', 60, ...
       'horizontalAlignment', 'center', 'verticalAlignment', 'middle', ...
       'fontSize', 16, 'fontWeight', 'bold')


% adjust_these lines
% $$$ m_line([lon_min, -16.77],[lat_max, 56.65],'color','r','linewi',1,'linestyle', '-');        
% $$$ m_line([lon_min, -16.77],[lat_min, 54.12],'color','r','linewi',1,'linestyle', '-');  
m_line([lon_max, -15.2],[lat_min, 57.2],'color','r','linewi',1,'linestyle', '-');        
m_line([lon_min, -19.7],[lat_min, 57.2],'color','r','linewi',1,'linestyle', '-');  
m_grid('box','fancy')
m_text(-14.5, 54.25, ['a) ' sprintf(datestr(time1))],'fontSize', 16, 'fontWeight', 'bold','BackgroundColor',[1 1 1]);

adjust_space


%% get bathym transect
m_line([-16.1120, -15.3720],[55.820, 55.023],'color','m','linewi',2,'linestyle', '--'); 
% $$$ [range,ln,lt]=m_lldist([-16.1120, -15.3720],[55.820, 55.023], 5000); %'m'
m_line([-15.969, -15.126],[55.820, 55.023],'color','g','linewi',2,'linestyle', '--'); 

%% Axis system
theta = -62;
U = [1 0];
V = [0 1];
[u,v] = rotate_vecd(U,V,theta);
m_vec(2,-14,55, u,v)
m_text(-13.5,55.3, 'y', 'fontSize', FS, 'fontWeight', 'bold')
m_text(-14.,54.75, 'x', 'fontSize', FS, 'fontWeight', 'bold')


%% --------- Large Area map (inset) -------- %
% $$$ lon_min=-32;
% $$$ lon_max=8;
% $$$ lat_min=42;
% $$$ lat_max=65;
lon_min=-32;
lon_max=8;
lat_min=42;
lat_max=65;

I=find(lat1m<lat_max & lat1m>lat_min);
J=find(lon1m<lon_max & lon1m>lon_min);
latitude=lat1m(I);
longitude=lon1m(J);
bathy=z1m(I,J);

%a2 = axes('position',[0.01 0.53 0.3 0.3]) ; % inset
a2 = axes('position',[0.03 0.17 0.2 0.2]) ; % inset
m_proj('mercator', 'long',[lon_min lon_max],'lat',[lat_min lat_max]);
[HH, HH] = m_contour(longitude,latitude,abs(bathy), [0:1000:4000], 'color', 'k');
m_gshhs_i('patch',[1 .9333 .6667]); %coastlines (Beige)                               
m_grid('box','fancy', 'yticklabels',[], 'xticklabels',[])

% Rectangle
lon_min=-20;
lon_max=-11;
lat_min=54;
lat_max=59;
[lomi, lami] = m_ll2xy(lon_min,lat_min);
[loma, lama] = m_ll2xy(lon_max,lat_max);
rectangle('Position', [lomi lami loma-lomi lama-lami],'linewidth', 1, 'edgecolor', 'r') 


%% Study region Map (zoom ++)
% $$$ lon_min=-15.85;
% $$$ lon_max=-15.75;
% $$$ lat_min=55.45;
% $$$ lat_max=55.5333;
lon_min=-15.85;
lon_max=-15.6;
lat_min=55.45;
lat_max=55.55;

I=find(lat<lat_max & lat>lat_min);
J=find(lon<lon_max & lon>lon_min);
latitude=lat(I);
longitude=lon(J);
bathy=z(I,J);

[lonLowRes, latLowRes] = meshgrid(longitude, latitude);
lonLowRes = lonLowRes(:);
latLowRes = latLowRes(:);
zLowRes = bathy(:);
I = find(zLowRes<-1000);
lonLowRes = lonLowRes(I);
latLowRes = latLowRes(I);
zLowRes = zLowRes(I);

% try to use high-res bathy (local refinement of gebco)
Y = lat_min:.001:lat_max;
X = lon_min:.001:lon_max;
[XI, YI] = meshgrid(X, Y);
%[XI, YI, ZI] = griddata(lonCWC, latCWC, zCWC, XI, YI);
[XI, YI, ZI] = griddata([lonCWC; lonLowRes], [latCWC; latLowRes], [zCWC; zLowRes], XI, YI);
latitude = YI;
longitude = XI;
bathy = ZI;


%a2 = axes('position',[0.01 0.16 0.25 0.35]) ; % inset
a3 = axes('position',[0.05 0.59 0.25 0.25]) ; % inset
m_proj('mercator', 'long',[lon_min lon_max],'lat',[lat_min lat_max]);
m_pcolor(longitude,latitude,abs(bathy)); shading flat;
V=[500:100:1200];
[HH, HH] = m_contour(longitude,latitude,abs(bathy), V, 'color', 'k');
caxis([min(V) max(V)])
m_grid('box','fancy', 'yticklabels',[], 'xticklabels',[])
set(HH, 'color', 'k', 'ShowText','on')%, 'LabelSpacing', 800)

% CTD casts
ctdLat  = load('/home/cyrf0006/research/NIOZ/RockallBank/CTD360/fred_processed/lat24h.txt');
ctdLon  = load('/home/cyrf0006/research/NIOZ/RockallBank/CTD360/fred_processed/lon24h.txt');
latVec = ctdLat(:,1) + ctdLat(:,2)/60;
lonVec = [ctdLon(:,1) + ctdLon(:,2)/60]*-1;
for i = 1:7%length(latVec)
    m_line(lonVec(i), latVec(i), 'marker', '.', 'color', 'm','markersize',14);
end

latVec2 = [55.5318; 55.5284; 55.5251; 55.5212; 55.5179; 55.5143];
lonVec2 = [-15.6643; -15.6604; -15.6565; -15.6527; -15.6491; -15.6454];
for i = 1:length(latVec2)
    m_line(lonVec2(i), latVec2(i), 'marker', '.', 'color', 'g','markersize',14);
end

% $$$ theta = -62;
% $$$ U = [1 0];
% $$$ V = [0 1];
% $$$ [u,v] = rotate_vecd(U,V,theta);
% $$$ m_vec(2,-15.79,55.47, u,v)

% $$$ m_text(-15.765,55.476, 'y', 'fontSize', FS, 'fontWeight', 'bold')
% $$$ m_text(-15.787,55.458, 'x', 'fontSize', FS, 'fontWeight', 'bold')
m_line([-16.1120, -15.3720],[55.820, 55.023],'color','m','linewi',1,'linestyle', '--'); 
m_line([-15.969, -15.126],[55.820, 55.023],'color','g','linewi',1,'linestyle', '--'); 

%get(a3, 'pos')


%% SUBPLOT 2
s = subplot(122);
[Y, I] = min(abs(timeVecTide-t0));
time2 = timeVecTide(I);
vMat = squeeze(tide.vv(:,:,I));
uMat = squeeze(tide.uu(:,:,I));
[uMat2, vMat2] = rotate_vecd(uMat, vMat, 62);
vMat = uMat2;

% refresh min/max
lon_min=-20;
lon_max=-11;
lat_min=54;
lat_max=59;

I=find(lat<lat_max & lat>lat_min);
J=find(lon<lon_max & lon>lon_min);
latitude=lat(I);
longitude=lon(J);
bathy=z(I,J);

m_proj('mercator','long',[min(longitude) max(longitude)],'lat',[min(latitude) max(latitude)]);

V=[0:100:4000];
m_contourf(lonVecTide,latVecTide,vMat/100, 50, 'lineStyle', 'none');
hold on
[HH, HH] = m_contour(longitude,latitude,abs(bathy), [0:500:2500], 'color', 'k');
m_gshhs_h('patch',[1 .9333 .6667]); %coastlines (Beige)                               
xlabel('Longitude', 'FontSize', FS, 'fontweight', 'bold')
set(gca, 'fontsize', FS)
caxis([-.25 .25])

% Rectangle around study region (zoom++ )
lon_min=-15.85;
lon_max=-15.6;
lat_min=55.45;
lat_max=55.55;

[lomi, lami] = m_ll2xy(lon_min,lat_min);
[loma, lama] = m_ll2xy(lon_max,lat_max);
rectangle('Position', [lomi lami loma-lomi lama-lami],'linewidth', 1, 'edgecolor', 'r') 
m_grid('box','fancy', 'yticklabel', [])
m_text(-14.5, 54.25, ['b) ' sprintf(datestr(time2))],'fontSize', 16, 'fontWeight', 'bold','BackgroundColor',[1 1 1]);

adjust_space



%% COLORBAR
c = colorbar('location','southoutside');
cb_pos = get(c, 'position');
cb_back = cb_pos;
cb_pos(4) = .02;
cb_pos(1) = .25;%cb_pos(1)+.1*cb_pos(3);
cb_pos(3) = .5;
cb_pos(2) = .9;%cb_pos(2)+cbar_offset;
set(c, 'pos', cb_pos);
set(c, 'fontsize', FS, 'fontweight', 'bold','color',[1 1 1])

ti = ylabel(c,'u (m s^{-1})', 'FontSize', FS, 'fontweight', 'bold','color',[1 1 1]*0);
ti_pos = get(ti, 'position');
set(ti, 'Rotation',0.0);     
clim = get(gca, 'clim');
set(ti, 'position', [0  2.5 ti_pos(3)]); 
set(gca, 'fontSize', FS)


%% SAVE FIGURE
set(gcf, 'renderer', 'painters'); % vectorial figure
print('-depsc', 'RockallMapAviso_subplot.eps')
