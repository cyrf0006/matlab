clear all
close all
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 2; % no. subplot column
nrow = 2; % no. subplot row
dx = 0.005 ; % horiz. space between subplots
dy = 0.01; % vert. space between subplots
lefs = 0.1; % very left of figure
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
paperwidth = 19;%cm
paperheight = 20;%cm
FS = 12;


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

% isolate 4 snaphots (to be edited)
[Y, I1] = min(abs(timeVecTide-datenum(2012,10,16,18,0,0)));
[Y, I2] = min(abs(timeVecTide-datenum(2012,10,17,0,0,0)));
[Y, I3] = min(abs(timeVecTide-datenum(2012,10,17,06,0,0)));
[Y, I4] = min(abs(timeVecTide-datenum(2012,10,17,12,0,0)));
vMat = squeeze(tide.vv(:,:,I1));
uMat = squeeze(tide.uu(:,:,I1));
[uMat1, vMat1] = rotate_vecd(uMat, vMat, 62);
vMat = squeeze(tide.vv(:,:,I2));
uMat = squeeze(tide.uu(:,:,I2));
[uMat2, vMat2] = rotate_vecd(uMat, vMat, 62);
vMat = squeeze(tide.vv(:,:,I3));
uMat = squeeze(tide.uu(:,:,I3));
[uMat3, vMat3] = rotate_vecd(uMat, vMat, 62);
vMat = squeeze(tide.vv(:,:,I4));
uMat = squeeze(tide.uu(:,:,I4));
[uMat4, vMat4] = rotate_vecd(uMat, vMat, 62);


%% --- HERE IS THE NEW CALCULATION --- %%
etaMat1 = squeeze(tide.slev(:,:,I1));
etaMat2 = squeeze(tide.slev(:,:,I2));
etaMat3 = squeeze(tide.slev(:,:,I3));
etaMat4 = squeeze(tide.slev(:,:,I4));

anom1 = nan(size(etaMat1));
anom2 = nan(size(etaMat1));
anom3 = nan(size(etaMat1));
anom4 = nan(size(etaMat1));
uanom1 = nan(size(etaMat1));
uanom2 = nan(size(etaMat1));
uanom3 = nan(size(etaMat1));
uanom4 = nan(size(etaMat1));
vanom1 = nan(size(etaMat1));
vanom2 = nan(size(etaMat1));
vanom3 = nan(size(etaMat1));
vanom4 = nan(size(etaMat1));

for i = 9:size(etaMat1,1)-9
    for j = 9:size(etaMat1,2)-9
        anom1(i,j) = etaMat1(i,j) - nanmean(nanmean(etaMat1(i-8:i+8, j-8:j+8)));
        anom2(i,j) = etaMat2(i,j) - nanmean(nanmean(etaMat2(i-8:i+8, j-8:j+8)));
        anom3(i,j) = etaMat3(i,j) - nanmean(nanmean(etaMat3(i-8:i+8, j-8:j+8)));
        anom4(i,j) = etaMat4(i,j) - nanmean(nanmean(etaMat4(i-8:i+8, j-8:j+8)));
        uanom1(i,j) = uMat1(i,j) - nanmean(nanmean(uMat1(i-8:i+8, j-8:j+8)));
        uanom2(i,j) = uMat2(i,j) - nanmean(nanmean(uMat2(i-8:i+8, j-8:j+8)));
        uanom3(i,j) = uMat3(i,j) - nanmean(nanmean(uMat3(i-8:i+8, j-8:j+8)));
        uanom4(i,j) = uMat4(i,j) - nanmean(nanmean(uMat4(i-8:i+8, j-8:j+8)));    
        vanom1(i,j) = vMat1(i,j) - nanmean(nanmean(vMat1(i-8:i+8, j-8:j+8)));
        vanom2(i,j) = vMat2(i,j) - nanmean(nanmean(vMat2(i-8:i+8, j-8:j+8)));
        vanom3(i,j) = vMat3(i,j) - nanmean(nanmean(vMat3(i-8:i+8, j-8:j+8)));
        vanom4(i,j) = vMat4(i,j) - nanmean(nanmean(vMat4(i-8:i+8, j-8:j+8)));        
    end
end
% ------------------------------------- %


%% ------- Plot main figure ------ %%
%h = figure('visible', 'off');
h = figure(1);
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 paperwidth paperheight])

load BR_symetric % <---- for velocities
colormap(BR_symetric);
%load ~/PhD/bathym/gebco64  <---- for bathymetry
%colormap(gebco);

for i = 1:4
    
% $$$     field = sprintf('uMat%d',i); %<-- here 'field' is a string
% $$$     uField = eval(field);
% $$$     field = sprintf('vMat%d',i); %<-- here 'field' is a string
% $$$     vField = eval(field);
    field = sprintf('uanom%d',i); %<-- here 'field' is a string
    uField = eval(field);
    field = sprintf('vanom%d',i); %<-- here 'field' is a string
    vField = eval(field);    
    field = sprintf('anom%d',i); %<-- here 'field' is a string
    aField = eval(field); %<-- now 'field' is the matrix!
    
    figure(1)
    s = subplot(2,2,i);
    m_proj('mercator','long',[min(longitude) max(longitude)],'lat',[min(latitude) max(latitude)]);
    %m_pcolor(longitude,latitude,abs(bathy)); shading flat;
    %shading flat
    m_contourf(lonVecTide,latVecTide,aField, 50, 'lineStyle', 'none');
    hold on
    [HH, HH] = m_contour(longitude,latitude,abs(bathy), [0:500:2500], 'color', 'k');
    
    %% Current arrows
    quiverDec = 5;
    Z = latVecTide(1:quiverDec:end);
    X = lonVecTide(1:quiverDec:end);
    U = uField(1:quiverDec:end, 1:quiverDec:end);
    V = vField(1:quiverDec:end, 1:quiverDec:end);
    [XX, ZZ] = meshgrid(X, Z);
    
    [kmax imax] = size(U);
    vecColor= [1 1 1]*0.5;
    
    m_quiver(XX, ZZ, U, V,6, 'color', 'k')
    
    
% $$$     %% 2D vector field
% $$$     scale = .3;
% $$$     for i = 1:imax
% $$$         for k = 1:kmax-1
% $$$             m_arrow(X(i),Z(k),U(k,i),V(k,i),scale);
% $$$         end
% $$$     end 
    caxis([-5 5])
    

    
    if i == 4 % < ---- For the colorbar
        c = colorbar('location','southoutside');
        cb_pos = get(c, 'position');
        cb_back = cb_pos;
        cb_pos(1) = .75; % x0
        cb_pos(2) = .15; % y0
        cb_pos(3) = .15; % width
        cb_pos(4) = .02; % thickness
        set(c, 'pos', cb_pos);
        set(c, 'fontsize', FS, 'fontweight', 'bold','color',[1 1 1])

        ti = ylabel(c,'\eta'' (cm)', 'FontSize', FS, 'fontweight', 'bold','color',[1 1 1]);
        ti_pos = get(ti, 'position');
        set(ti, 'Rotation',0.0);     
        clim = get(gca, 'clim');
        set(ti, 'position', [0  2.5 ti_pos(3)]); 
        set(s, 'fontSize', FS)
    end
    
    adjust_space
    
    if i == 1 
        m_grid('box','fancy', 'xticklabels',[])
        ylabel('Longitude', 'FontSize', FS, 'fontweight', 'bold')
    elseif i == 2
        m_grid('box','fancy', 'yticklabels',[], 'xticklabels',[])
    elseif i == 3
        m_grid('box','fancy')
        ylabel('Longitude', 'FontSize', FS, 'fontweight', 'bold')
        xlabel('Latitude', 'FontSize', FS, 'fontweight', 'bold')
    elseif i == 4
        m_grid('box','fancy', 'yticklabels',[])
        xlabel('Latitude', 'FontSize', FS, 'fontweight', 'bold')
    end
    

    lon_min2=-15.85;
    lon_max2=-15.75;
    lat_min2=55.45;
    lat_max2=55.5333;

    [lomi, lami] = m_ll2xy(lon_min2,lat_min2);
    [loma, lama] = m_ll2xy(lon_max2,lat_max2);
    rectangle('Position', [lomi lami loma-lomi lama-lami],'linewidth', 1, 'edgecolor', 'r') 
    
    if i == 1  % <----- Big 'if' statement that you can remove...
        % Rectangle around study region (zoom++ )
% $$$         lon_min=-15.85;
% $$$         lon_max=-15.75;
% $$$         lat_min=55.45;
% $$$         lat_max=55.5333;
% $$$ 
% $$$         [lomi, lami] = m_ll2xy(lon_min,lat_min);
% $$$         [loma, lama] = m_ll2xy(lon_max,lat_max);
% $$$         rectangle('Position', [lomi lami loma-lomi lama-lami],'linewidth', 1, 'edgecolor', 'r') 

        % few annotations
        m_text(-14.75, 56.6667, 'Rockall Bank', 'rotation', 60, ...
               'horizontalAlignment', 'center', 'verticalAlignment', 'middle', ...
               'fontSize', 14, 'fontWeight', 'bold')


% $$$         % adjust_these lines
% $$$         m_line([lon_min, -14.1],[lat_max, 56.55],'color','r','linewi',1,'linestyle', '-');        
% $$$         m_line([lon_min, -14.1],[lat_min, 54.15],'color','r','linewi',1,'linestyle', '-');  
% $$$         %        m_grid('box','fancy')
% $$$ 
% $$$ 
% $$$         %% get bathym transect
% $$$         m_line([-16.1120, -15.3720],[55.820, 55.023],'color','m','linewi',2,'linestyle', '--'); 
% $$$         [range,ln,lt]=m_lldist([-16.1120, -15.3720],[55.820, 55.023], 5000); %'m'

        % put the axis system
        theta = -62;
        U = [1 0];
        V = [0 1];
        [u,v] = rotate_vecd(U,V,theta);
        %m_vec(2,-15.775,55.455, u,v)
        m_vec(2,-13,55, u,v)
        m_text(-12,55.6, 'y', 'fontSize', FS, 'fontWeight', 'bold')
        m_text(-13,54.5, 'x', 'fontSize', FS, 'fontWeight', 'bold')
        
        %% --------- Large Area map (inset) -------- %
        lon_min=-32;
        lon_max=8;
        lat_min=42;
        lat_max=65;

        I=find(lat1m<lat_max & lat1m>lat_min);
        J=find(lon1m<lon_max & lon1m>lon_min);
        latitudeLarge=lat1m(I);
        longitudeLarge=lon1m(J);
        bathyLarge = z1m(I,J);

        a2 = axes('position',[0.13 0.78 0.15 0.15]) ; % inset
        m_proj('mercator', 'long',[lon_min lon_max],'lat',[lat_min lat_max]);
        [HH, HH] = m_contour(longitudeLarge,latitudeLarge,abs(bathyLarge), [0:1000:4000], 'color', 'k');
        m_gshhs_c('patch',[1 .9333 .6667]); %coastlines (Beige)                               
        m_grid('box','fancy', 'yticklabels',[], 'xticklabels',[])

        % Rectangle
        lon_min=-20;
        lon_max=-11;
        lat_min=54;
        lat_max=59;
        [lomi, lami] = m_ll2xy(lon_min,lat_min);
        [loma, lama] = m_ll2xy(lon_max,lat_max);
        rectangle('Position', [lomi lami loma-lomi lama-lami],'linewidth', 1, 'edgecolor', 'r') 


% $$$         %% Study region Map (zoom ++)
% $$$         lon_min=-15.85;
% $$$         lon_max=-15.75;
% $$$         lat_min=55.45;
% $$$         lat_max=55.5333;
% $$$ 
% $$$         I=find(lat<lat_max & lat>lat_min);
% $$$         J=find(lon<lon_max & lon>lon_min);
% $$$         latitudeZoom=lat(I);
% $$$         longitudeZoom=lon(J);
% $$$         bathyZoom=z(I,J);
% $$$ 
% $$$         a2 = axes('position',[0.31 0.54 0.2 0.2]) ; % inset
% $$$         m_proj('mercator', 'long',[lon_min lon_max],'lat',[lat_min lat_max]);
% $$$         %m_contourf(longitude,latitude,abs(bathy),V, 'linestyle', 'none');
% $$$         %        m_pcolor(longitudeZoom,latitudeZoom,abs(bathyZoom)); shading flat;
% $$$         V=[500:100:1200];
% $$$         [HH, HH] = m_contour(longitudeZoom,latitudeZoom,abs(bathyZoom), V, 'color', 'k');
% $$$         caxis([min(V) max(V)])
% $$$         m_grid('box','fancy', 'yticklabels',[], 'xticklabels',[])
% $$$         set(HH, 'color', 'k', 'ShowText','on')%, 'LabelSpacing', 800)
% $$$ 
% $$$ 
% $$$         % Example to put CTD casts location
% $$$ % $$$         % CTD casts
% $$$ % $$$         ctdLat  = load('/home/cyrf0006/research/NIOZ/RockallBank/CTD360/fred_processed/lat24h.txt');
% $$$ % $$$         ctdLon  = load('/home/cyrf0006/research/NIOZ/RockallBank/CTD360/fred_processed/lon24h.txt');
% $$$ % $$$         latVec = ctdLat(:,1) + ctdLat(:,2)/60;
% $$$ % $$$         lonVec = [ctdLon(:,1) + ctdLon(:,2)/60]*-1;
% $$$ % $$$         for i = 1:7%length(latVec)
% $$$ % $$$             m_line(lonVec(i), latVec(i), 'marker', '.', 'color', 'm','markersize',14);
% $$$ % $$$         end
% $$$ 
% $$$         % put the axis system
% $$$         theta = -62;
% $$$         U = [1 0];
% $$$         V = [0 1];
% $$$         [u,v] = rotate_vecd(U,V,theta);
% $$$         %m_vec(2,-15.775,55.455, u,v)
% $$$         m_vec(2,-15.79,55.47, u,v)
% $$$ 
% $$$         % Mooring location
% $$$         m_line(-15.7975, 55.4824,'marker','p','color',[0 0 0], 'markerfacecolor', 'r','markersize',14);
% $$$         %m_line(-15.8, 55.49,'marker','p','color',[0 0 0], 'markerfacecolor', 'm','markersize',10);
% $$$         m_text(-15.765,55.476, 'y', 'fontSize', FS, 'fontWeight', 'bold')
% $$$         m_text(-15.787,55.458, 'x', 'fontSize', FS, 'fontWeight', 'bold')
% $$$         m_line([-16.1120, -15.3720],[55.820, 55.023],'color','m','linewi',1,'linestyle', '--'); 
    end
    
    
    
end

% Save figure in eps (renderer 'painters' means vectorial)
set(gcf, 'renderer', 'painters'); % vectorial figure
print('-depsc', 'AvisoMapExample.eps')
