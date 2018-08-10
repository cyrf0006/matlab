clear all
close all


%% plotting info
%figure dimension
paperwidth = 20;%cm
paperheight = 30;%cm
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
% $$$ % Lower res.
% $$$ %decimFactor = 4;
% $$$ decimFactor = 40;
% $$$ lat1m = lat(1:decimFactor:end);
% $$$ lon1m = lon(1:decimFactor:end);
% $$$ z1m = z(1:decimFactor:end, 1:decimFactor:end);

%% Get local bathym (multibeam)
path2 = ['~/research/NIOZ/RockallBank/bathymetry/XYZHERMESALL_ROCKALL.CSV'];
lims = [40 70 -35 10];
CWC = load(path2);

latCWC = CWC(:,1);
lonCWC = CWC(:,2);
zCWC = CWC(:,3);
%reduce size (zoom ++)
% $$$ lon_min=-15.85;
% $$$ lon_max=-15.6;
% $$$ lat_min=55.45;
% $$$ lat_max=55.55;
lon_min=min(lonCWC);
lon_max=max(lonCWC);
lat_min=min(latCWC);
lat_max=max(latCWC);
%lat_max=55.53;


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

I=find(latCWC<=lat_max & latCWC>=lat_min & lonCWC<=lon_max & lonCWC>=lon_min);
latCWC = latCWC(I);
lonCWC = lonCWC(I);
zCWC = zCWC(I);





%% ------- Plot main figure ------ %%
%h = figure('visible', 'off');
h = figure(1);
clf
%set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 paperwidth paperheight])

% try to use high-res bathy (local refinement of gebco)
Y = lat_min:.001:lat_max;
X = lon_min:.001:lon_max;
Y = lat_min:.05:lat_max;
X = lon_min:.05:lon_max;
[XI, YI] = meshgrid(X, Y);
%[XI, YI, ZI] = griddata(lonCWC, latCWC, zCWC, XI, YI);
[XI, YI, ZI] = griddata([lonCWC; lonLowRes], [latCWC; latLowRes], [zCWC; zLowRes], XI, YI);
latitude = YI;
longitude = XI;
bathy = ZI;

surf(longitude,latitude,abs(bathy), 'linestyle', 'none');

colormap(flipud(jet))

%V=[500:100:1200];
%[HH, HH] = m_contour(longitude,latitude,abs(bathy), V, 'color', 'k');
%caxis([min(V) max(V)])
%m_grid('box','fancy', 'yticklabels',[], 'xticklabels',[])
%set(HH, 'color', 'k', 'ShowText','on')%, 'LabelSpacing', 800)
set(gca, 'zdir', 'reverse')
set(gca, 'tickdir', 'in')
%set(gca, 'box', 'on')
% $$$ set(gca, 'yticklabel', [])
% $$$ set(gca, 'xticklabel', [])

xlim([lon_min lon_max])
ylim([lat_min lat_max])
zlim([500 2000])

hold on
% CTD casts
ctdLat  = load('/home/cyrf0006/research/NIOZ/RockallBank/CTD360/fred_processed/lat24h.txt');
ctdLon  = load('/home/cyrf0006/research/NIOZ/RockallBank/CTD360/fred_processed/lon24h.txt');
latVec = ctdLat(:,1) + ctdLat(:,2)/60;
lonVec = [ctdLon(:,1) + ctdLon(:,2)/60]*-1;
% $$$ for i = 1:7%length(latVec)
% $$$     [Y, I] = min(abs(longitude(:)-lonVec(i)) + abs(latitude(:)-latVec(i)));
% $$$     plot3(lonVec(i), latVec(i), abs(bathy(I)),'marker', 'o', 'color', ...
% $$$           [1 1 1],'markersize',14);
% $$$ end
% $$$ 
% $$$ latVec2 = [55.5318; 55.5284; 55.5251; 55.5212; 55.5179; 55.5143];
% $$$ lonVec2 = [-15.6643; -15.6604; -15.6565; -15.6527; -15.6491; -15.6454];
% $$$ for i = 1:length(latVec2)
% $$$     [Y, I] = min(abs(longitude(:)-lonVec2(i)) + abs(latitude(:)-latVec2(i)));
% $$$     plot3(lonVec2(i), latVec2(i), abs(bathy(I)),'marker', 'o', ...
% $$$           'color', [1 1 1],'markersize',14);
% $$$ end


[range,ln,lt]=m_lldist([-16.1120, -15.3720],[55.820, 55.023], 1000); %'m'
zn = nan(size(lt));
for i = 1:length(lt)
    [Y, I] = min(abs(longitude(:)-ln(i)) + abs(latitude(:)-lt(i)));
    zn(i) = abs(bathy(I));
end
plot3(ln, lt, zn, 'color','m','linewi',3,'linestyle', '-')

[range,ln,lt]=m_lldist([-15.969, -15.126],[55.820, 55.023], 1000); %'m'
zn = nan(size(lt));
for i = 1:length(lt)
    [Y, I] = min(abs(longitude(:)-ln(i)) + abs(latitude(:)-lt(i)));
    zn(i) = abs(bathy(I));
end
plot3(ln, lt, zn, 'color','g','linewi',3,'linestyle', '-')

whitebg([0 0 0])  

% $$$ set(gca, 'cameraPosition', [-0.0157 0.0553 -5.0321]*1000)
% $$$ set(gca, 'cameraViewAngle', 8.4)
% $$$ set(gca, 'cameraTarget', [-15.725 55.5 850])


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

ti = ylabel(c,'Depth (m)', 'FontSize', FS, 'fontweight', 'bold','color',[1 1 1]*0);
ti_pos = get(ti, 'position');
set(ti, 'Rotation',0.0);     
clim = get(gca, 'clim');
set(ti, 'position', [0  ti_pos(2) ti_pos(3)]); 
set(gca, 'fontSize', FS)


%% SAVE FIGURE
set(gcf, 'renderer', 'painters'); % vectorial figure
print('-depsc', 'Rockall_multiBeam.eps')
