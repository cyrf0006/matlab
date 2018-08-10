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
paperwidth = 22;%cm
paperheight = 20;%cm
FS = 12;

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


%% Get local bathym (multibeam)
path2 = ['~/research/NIOZ/RockallBank/bathymetry/coral_XYZ_clean.csv'];
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
lon_max=-10;
lat_min=54;
lat_max=60;

I=find(lat<lat_max & lat>lat_min);
J=find(lon<lon_max & lon>lon_min);
latitude=lat(I);
longitude=lon(J);
bathy=z(I,J);


%% Plot map
%h = figure('visible', 'off');
h = figure(1);
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 paperwidth paperheight])
m_proj('mercator','long',[min(longitude) max(longitude)],'lat',[min(latitude) max(latitude)]);
load ~/PhD/bathym/gebco64
colormap(gebco);

V=[0:100:4000];
m_pcolor(longitude,latitude,abs(bathy)); shading flat;
hold on
[HH, HH] = m_contour(longitude,latitude,abs(bathy), [0:500:4000], 'color', 'k');
m_gshhs_c('patch',[1 .9333 .6667]); %coastlines (Beige)                               
xlabel('Longitude', 'FontSize', FS, 'fontweight', 'bold')
ylabel('Latitude', 'FontSize', FS, 'fontweight', 'bold')
caxis([min(V) max(V)])
%set(gca, 'fontsize', FS)

% Rectangle around study region (zoom++ )
lon_min=-15.85;
lon_max=-15.6;
lat_min=55.45;
lat_max=55.55;

[lomi, lami] = m_ll2xy(lon_min,lat_min);
[loma, lama] = m_ll2xy(lon_max,lat_max);
rectangle('Position', [lomi lami loma-lomi lama-lami],'linewidth', 1, 'edgecolor', 'r') 

% few annotations
m_text(-14, 57.2, 'Rockall Bank', 'rotation', 66, ...
       'horizontalAlignment', 'center', 'verticalAlignment', 'middle', ...
       'fontSize', 14, 'fontWeight', 'bold')


% adjust_these lines
% $$$ m_line([lon_min, -13.4],[lat_max, 55.95],'color','r','linewi',1,'linestyle', '-');        
% $$$ m_line([lon_min, -13.4],[lat_min, 54.15],'color','r','linewi',1,'linestyle', '-');  
m_line([lon_max, -13.88],[lat_min, 57.8],'color','r','linewi',1,'linestyle', '-');        
m_line([lon_min, -19.58],[lat_min, 57.8],'color','r','linewi',1,'linestyle', '-');  
m_grid('box','fancy')


%% get bathym transect

% --------> T1 (Haas Mound) <--------- %
m_line([-16.1120, -15.3720],[55.820, 55.023],'color','m','linewi',2,'linestyle', '--'); 
[range,ln,lt]=m_lldist([-16.1120, -15.3720],[55.820, 55.023], 5000); %'m'

% distance_vector
dv = [0:range./length(ln):range-range./length(ln)]';
% depth vector
zv = nan(length(dv),1);
% depth at transect points
I = find(lon>=min(ln) & lon<=max(max(ln)));
J = find(lat>=min(lt) & lat<=max(max(lt)));
[X, Y] = meshgrid(lon(I), lat(J));
All_pts = [X(:) Y(:)]; %very all pts of the map
bathy = z(J,I);
bathy = bathy(:);
clear X Y

% **low resolution**
for i=1:length(ln)      
    % find least square (closest point)
    [Y, I] = min((ln(i)-All_pts(:,1)).^2+(lt(i)-All_pts(:,2)).^2);
    % depth vector corresponding to distance vector
    zv(i) = bathy(I);
end
transect = [ln' lt' dv zv];
dlmwrite('RockallBathym_CTD.dat', transect ,'delimiter',' ','precision',12) 

% **high resolution**
[rangeHR,lnHR,ltHR]=m_lldist([-15.85, -15.767],[55.539, 55.45], 5000); %'m'
% note that these pts are only used to compute depth, not to plot
% on map (otherwise 500pts are 'diluted' over a wider distance)
dvHR = [0:rangeHR./length(lnHR):rangeHR-rangeHR./length(lnHR)]';
zvHR = nan(length(dvHR),1);
dist = nan(length(dvHR),1);
for i=1:length(lnHR)      
    % find least square (closest point)
    [Y, I] = min(abs(lnHR(i)-lonCWC)+abs(ltHR(i)-latCWC));
    dist(i) = Y;
    % depth vector corresponding to distance vector
    zvHR(i) = zCWC(I);
end
shift = m_lldist([ln(1) lnHR(1)], [lt(1) ltHR(1)],1);
dvHR = dvHR+shift;
transect = [lnHR' ltHR' dvHR zvHR];
dlmwrite('RockallBathym_CWC.dat', transect ,'delimiter',' ','precision',12) 

% --------> T2  (smaller Mound) <-------- %
m_line([-15.969, -15.126],[55.820, 55.023],'color','g','linewi',2,'linestyle', '--'); 
[range,ln,lt]=m_lldist([-15.969, -15.126],[55.820, 55.023], 5000); 

% distance_vector
dv = [0:range./length(ln):range-range./length(ln)]';
% depth vector
zv = nan(length(dv),1);
% depth at transect points
I = find(lon>=min(ln) & lon<=max(max(ln)));
J = find(lat>=min(lt) & lat<=max(max(lt)));
[X, Y] = meshgrid(lon(I), lat(J));
All_pts = [X(:) Y(:)]; %very all pts of the map
bathy = z(J,I);
bathy = bathy(:);
clear X Y

% **low resolution**
for i=1:length(ln)      
    % find least square (closest point)
    [Y, I] = min((ln(i)-All_pts(:,1)).^2+(lt(i)-All_pts(:,2)).^2);
    % depth vector corresponding to distance vector
    zv(i) = bathy(I);
end
transect = [ln' lt' dv zv];
dlmwrite('RockallBathym_CTD_2.dat', transect ,'delimiter',' ','precision',12) 


% **high resolution**
[rangeHR,lnHR,ltHR]=m_lldist([-15.7375, -15.6],[55.6, 55.4750], 5000);%'g'
dvHR = [0:rangeHR./length(lnHR):rangeHR-rangeHR./length(lnHR)]';
zvHR = nan(length(lnHR),1);
dist = nan(length(lnHR),1);
for i=1:length(lnHR)      
    % find least square (closest point)
    [Y, I] = min(abs(lnHR(i)-lonCWC)+abs(ltHR(i)-latCWC));
    dist(i) = Y;
    % depth vector corresponding to distance vector
    zvHR(i) = zCWC(I);
end
% remove depth outside multibeam (visual inspection of 'dist')
zvHR(1:2000) = NaN;
zvHR(4260:end) = NaN;

shift = m_lldist([ln(1) lnHR(1)], [lt(1) ltHR(1)],1);
dvHR = dvHR+shift;
transect = [lnHR' ltHR' dvHR zvHR];
dlmwrite('RockallBathym_CWC_2.dat', transect ,'delimiter',' ','precision',12) 


%% Large Area map (inset)
lon_min=-32;
lon_max=8;
lat_min=42;
lat_max=65;

I=find(lat1m<lat_max & lat1m>lat_min);
J=find(lon1m<lon_max & lon1m>lon_min);
latitude=lat1m(I);
longitude=lon1m(J);
bathy=z1m(I,J);

a2 = axes('position',[0.57 0.13 0.28 0.28]) ; % inset
m_proj('mercator', 'long',[lon_min lon_max],'lat',[lat_min lat_max]);
hold on
m_pcolor(longitude,latitude,abs(bathy)); shading flat;
[HH, HH] = m_contour(longitude,latitude,abs(bathy), [0:1000:4000], 'color', 'k');
m_gshhs_i('patch',[1 .9333 .6667]); %coastlines (Beige)                               
caxis([min(V) max(V)])
m_grid('box','fancy', 'yticklabels',[], 'xticklabels',[])

load ~/PhD/bathym/gebco64
colormap(gebco);

% Rectangle
lon_min=-20;
lon_max=-10;
lat_min=54;
lat_max=60;
[lomi, lami] = m_ll2xy(lon_min,lat_min);
[loma, lama] = m_ll2xy(lon_max,lat_max);

rectangle('Position', [lomi lami loma-lomi lama-lami],'linewidth', 1, 'edgecolor', 'r') 

cbar_offset = 0;%0.04; 
offset2 =0;% 0.02; % offset between heigth of colorbar 
c = colorbar('location','southoutside');
cb_pos = get(c, 'position');
cb_back = cb_pos;
cb_pos(4) = .02;
cb_pos(1) = .2;%cb_pos(1)+.1*cb_pos(3);
cb_pos(3) = .25;
cb_pos(2) = .15;%cb_pos(2)+cbar_offset;
set(c, 'pos', cb_pos);
set(c, 'fontsize', FS, 'fontweight', 'bold')

ti = ylabel(c,'Depth (m)', 'FontSize', FS, 'fontweight', 'bold');
ti_pos = get(ti, 'position');
set(ti, 'Rotation',0.0);     
clim = get(gca, 'clim');
set(ti, 'position', [(clim(2)-clim(1))./2  2 ti_pos(3)]); 
set(gca, 'fontSize', FS)



%% Study region Map (zoom ++)
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

%a2 = axes('position',[0.55 0.12 0.4 0.3]) ; % inset
a2 = axes('position',[0.2 0.61 0.4 0.3]) ; % inset
m_proj('mercator', 'long',[lon_min lon_max],'lat',[lat_min lat_max]);
%m_contourf(longitude,latitude,abs(bathy),V, 'linestyle', 'none');
%m_pcolor(longitude,latitude,abs(bathy)); shading flat;
V=[600 850 1000 1200];
[HH, HH] = m_contour(longitude,latitude,abs(bathy), V, 'color', 'k');
%caxis([min(V) max(V)])
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
m_line([-15.969, -15.126],[55.820, 55.023],'color','g','linewi',1,'linestyle', '--'); 


theta = -62;
U = [1 0];
V = [0 1];
[u,v] = rotate_vecd(U,V,theta);
%m_vec(2,-15.775,55.455, u,v)
m_vec(2,-15.79,55.47, u,v)


m_line(-15.7975, 55.4824,'marker','p','color',[0 0 0], 'markerfacecolor', 'r','markersize',14);
%m_line(-15.8, 55.49,'marker','p','color',[0 0 0], 'markerfacecolor', 'm','markersize',10);
m_text(-15.765,55.476, 'y', 'fontSize', FS, 'fontWeight', 'bold')
m_text(-15.787,55.458, 'x', 'fontSize', FS, 'fontWeight', 'bold')
m_line([-16.1120, -15.3720],[55.820, 55.023],'color','m','linewi',1,'linestyle', '--'); 



%save figure
%print('-dpng', '-r300', 'RockallMap.png')
set(gcf, 'renderer', 'painters'); % vectorial figure
print('-depsc', 'RockallMap.eps')

