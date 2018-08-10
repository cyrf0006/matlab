function fr_chla_map(mapFile, varargin)    

%function fr_chla_map(mapFile, varargin)    
%
% To be run in /home/cyrf0006/IML/fronts/Chl-a
%
%
% usage ex:
% fr_chla_map('Labrador.hdf', 'East');
%   OR
% fr_chla_map('Labrador.hdf')
%
% Running without varargin will load previous display
% setup. Running with varargin will ERASE previous display setup.
%

% other ex;
% fr_chla_map('Hudson-sud.hdf', 'North');
% fr_chla_map('HudsonStrait.hdf', 'North');
% fr_chla_map('Hudson-est.hdf', 'East');   
% etc.
% 
% F. Cyr - August 2013

% You can also give the complete path so no need to be in a
% specific folder to run it ()
% note from F. Cyr - November 2014
% ex: fr_chla_map('~/research/fronts/Chl-a/Golfe.hdf');   

%
    
% Get corners (given by P. Larouche)
region = mapFile(1:end-4);    
I = find(region=='/');
if ~isempty(I)
    region(1:I(end)) = [];
end

if strcmp(region, 'Alaska')    
    lonMin = -150; latMax = 62; lonMax = -130; latMin = 50;
    projLon = [lonMin lonMax]; projLat = [latMin latMax]; % for projection Lims, adjust here
elseif strcmp(region, 'Baffin')
    lonMin = -84; latMax = 82; lonMax = -50; latMin = 60;
    projLon = [lonMin lonMax]; projLat = [latMin latMax];
elseif strcmp(region, 'Golfe')
    lonMin = -70; latMax = 52; lonMax = -55; latMin = 45;
    projLon = [lonMin lonMax]; projLat = [latMin latMax];
elseif strcmp(region, 'Hudson-est')
    lonMin = -84; latMax = 59; lonMax = -76; latMin = 51;    
    projLon = [lonMin lonMax]; projLat = [latMin latMax];
elseif strcmp(region, 'Hudson-ouest')
    lonMin = -95; latMax = 64; lonMax = -90; latMin = 60;
    projLon = [lonMin lonMax]; projLat = [latMin latMax];
elseif strcmp(region, 'Hudson-sud')
    lonMin = -96; latMax = 60; lonMax = -82; latMin = 55;    
    projLon = [lonMin lonMax]; projLat = [latMin latMax];
elseif strcmp(region, 'Hudson')
    lonMin = -96; latMax = 71; lonMax = -72; latMin = 51;    
    projLon = [lonMin lonMax]; projLat = [latMin latMax];
elseif strcmp(region, 'HudsonStrait')
    lonMin = -80; latMax = 65; lonMax = -60; latMin = 58;    
    projLon = [lonMin lonMax]; projLat = [latMin latMax];
elseif strcmp(region, 'Labrador')
    lonMin = -64; latMax = 59; lonMax = -59; latMin = 52;    
    projLon = [lonMin lonMax]; projLat = [latMin latMax];
elseif strcmp(region, 'NovaScotia')
    lonMin = -72; latMax = 47; lonMax = -57; latMin = 40;    
    projLon = [lonMin -58]; projLat = [latMin latMax];
elseif strcmp(region, 'NFDL')
    lonMin = -60; latMax = 52; lonMax = -44; latMin = 44;    
    %    projLon = [lonMin -51]; projLat = [latMin latMax];
    projLon = [lonMin lonMax]; projLat = [latMin latMax];
elseif strcmp(region, 'Vancouver')
    lonMin = -132; latMax = 55; lonMax = -122; latMin = 44;        
    projLon = [lonMin lonMax]; projLat = [latMin latMax];
elseif strcmp(region, 'Pacifique')
    lonMin = -160; latMax = 62; lonMax = -120; latMin = 40;        
    projLon = [-157 -126]; projLat = [50 latMax];
else
    disp('Have wrong corner lat-lon!')
    return
end

% Get data and build lat lon vectors
data = hdfread(mapFile,'l3m_data');
I = find(data==-32767);
data(I) = NaN;
data = double(data);
dlat = (latMax-latMin)/(size(data,1)-1);
dlon = (lonMax-lonMin)/(size(data,2)-1);
latVec = latMin:dlat:latMax;
lonVec = lonMin:dlon:lonMax;

% Check for figure dimensions, colorbar position, etc.
% *********************** Adjust_space.m ************************ %
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
count_col = 1;
count_row = 1;
% *************************************************************** %
if ~isempty(varargin) == 1 
    if strcmp(varargin{1}, 'North') == 1
        paperwidth = 23;%cm
        paperheight = 17;%cm
        cbar_width = 0.02;
        cbar_offset = 0.16; % colorbar offset from figure
        offset2 = 0.1; % offset between heigth of colorbar 
        tipos1 = 0.15;
        lefs = 0.12; % very left of figure
        rigs = 0.05; % very right of figure
        tops = 0.1; % top of figure
        bots = 0.07; % bottom of figure
        figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
        figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
        cbar_loc = 'NorthOutside';
        north = 1;
    else        
        % (East)
        paperwidth = 17;%cm
        paperheight = 30;%cm
        cbar_width = 0.02;
        cbar_offset = 0.12; % colorbar offset from figure
        offset2 = 0.1; % offset between heigth of colorbar 
        tipos1 = 0.15;
        cbar_loc = 'EastOutside';
        lefs = 0.1; % very left of figure
        rigs = 0.12; % very right of figure
        tops = 0.05; % top of figure
        bots = 0.1; % bottom of figure
        figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
        figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
        north = 0;
    end    
else
    % If second varargin not provided, we assume these parameters
    % exists
    path_display = ['~/data/matlab_bathym/chla_maps/' region '_display.mat'];
    %path_display = ['~/data/matlab_bathym/' region '_display.mat'];
    load(path_display)
end


% Load bathymetry
path = ['~/data/matlab_bathym/chla_maps/' region '.mat'];
fid = fopen(path);
if fid == -1 % doesnt exist yet!
    disp('Region doesn''t exist yet!')
    disp(' -> Extract bathym from Gebco 30'''', may take some time...')
    lims = [latMin latMax lonMin lonMax];
    [lat lon z] = getGebco('~/data/GEBCO/GEBCO_08.nc', 'z', lims);
    disp('   done!')
    save(path, 'lat', 'lon', 'z');
else
    load(path)
end

% Now plot
h = figure('Visible', 'off');
%h=figure(1);
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 paperwidth paperheight])
m_proj('mercator','long',projLon,'lat', projLat);
%m_proj('mercator','long',[min(lonVec) max(lonVec)],'lat',[min(latVec) max(latVec)]);

m_grid('box','fancy')
hold on
%  WATCH OUT! make sure m_gshhs_h is  'high res' for final figure
m_pcolor(lonVec, latVec, log10(data)); shading interp
%m_gshhs_h('patch',[.7 .7 .7]); %coastlines    
m_gshhs_h('patch',[1 .9333 .6667]); %coastlines (Beige)                               
m_grid('box','fancy')

plot_border('~/research/fronts/borders/NorthAmerica.dat',[min(lon) max(lon)], [min(lat) max(lat)]);


% ---------- Xtra features on different plots ------------- %
if strcmp(region, 'Baffin') == 1
    [H3, H3] = m_contour(lon,lat,z, [-500 -1000], 'color', 'k');
    fr_putPlaces('BAFF')
    m_line(-69.3457, 77.4830, 'marker','.','MarkerFaceColor', ...
           'k','markersize',14,'color',[.6 0 0]);
end

if strcmp(region, 'NovaScotia') == 1
    [H3, H3] = m_contour(lon,lat,z, [-50 -100 -200 -1000], 'color', 'k');
    fr_putPlaces('SCOS')
end

if strcmp(region, 'NFDL') == 1 | strcmp(region, 'Labrador') == 1
    [H3, H3] = m_contour(lon,lat,z, [-100 -200 -1000], 'color', 'k');
    fr_putPlaces('NFDL')
end

if strcmp(region, 'Golfe') == 1
    [H3, H3] = m_contour(lon,lat,z, [-50 -300], 'color', 'k');
    fr_putPlaces('GULF')
    m_line(-69.716606, 48.156604,'marker','.','MarkerFaceColor', ...
           'k','markersize',14,'color',[.6 0 0]);
    m_line(-67.386367, 49.35, 'marker','.','MarkerFaceColor', ...
           'k','markersize',14,'color',[.6 0 0]);
end

if strcmp(region, 'Hudson') == 1 
    [H3, H3] = m_contour(lon,lat,z, [-100 -100], 'color', 'k');
    fr_putPlaces('HUDB')
end

if strcmp(region, 'HudsonStrait') == 1 
    [H3, H3] = m_contour(lon,lat,z, [-50 -200 -1000], 'color', 'k');
    fr_putPlaces('HUDS')
end

if strcmp(region, 'Hudson-est') == 1 | strcmp(region, 'Hudson-sud') == ...
        1 | strcmp(region, 'Hudson-ouest') == 1
    [H3, H3] = m_contour(lon,lat,z, [-25 -50 -100], 'color', 'k');
end

if strcmp(region, 'Alaska') == 1
    [H3, H3] = m_contour(lon,lat,z, [-100 -1000], 'color', 'k');
end

if strcmp(region, 'Pacifique') == 1 % in fact it is Gulf Alaska
    [H3, H3] = m_contour(lon,lat,z, [-100 -1000], 'color', 'k');
    fr_putPlaces('GALA')
end
 
if strcmp(region, 'Vancouver') == 1
    [H3, H3] = m_contour(lon,lat,z, [-100 -1000], 'color', 'k');
    fr_putPlaces('VANC')
end


% --------------------------------------------------------- %
hold off

xlabel('Longitude', 'FontSize', 10)
ylabel('Latitude', 'FontSize', 10)
set(gca, 'fontsize', 10)

adjust_space

%caxis([0 5])
if strcmp(region, 'NFDL') == 1 | strcmp(region, 'HudsonStrait') == ...
        1 | strcmp(region, 'Baffin') == 1 
    caxis([-1 1])
    YTICKS = [-1:.2:1];
else
    caxis([-1 1.2])
    YTICKS = [-1:.2:1.2];
end

if north
    c = colorbar(cbar_loc);
    Pos = get(gca, 'position');
    set(c, 'FontSize', 10, 'position', [Pos(1)+offset2 Pos(2)+Pos(4)+cbar_offset Pos(3)-2*offset2 cbar_width]);
    ti = ylabel(c,'log [chl-a (mg m^{-3})]', 'FontSize', 10);
    ti_pos = get(ti, 'position');
    set(ti, 'Rotation',0.0);     
    clim = get(gca, 'clim');
    set(ti, 'position', [(clim(2)-abs(clim(1)))./2  8 ti_pos(3)]); 
    set(c, 'xtick', YTICKS)
else
    c = colorbar(cbar_loc);
    Pos = get(gca, 'position');
    set(c, 'FontSize', 10, 'position', [Pos(1)+Pos(3)+cbar_offset Pos(2)+offset2 cbar_width Pos(4)-2*offset2]);
    ti = ylabel(c,'log [chl-a (mg m^{-3})]', 'FontSize', 10);
    ti_pos = get(ti, 'position');
    set(ti, 'position', [ti_pos(1)-tipos1 ti_pos(2) ti_pos(3)]);
    set(c, 'ytick', YTICKS)
end


%  WATCH OUT! Here -r100 is quite low resolution, just to speed up
%  treatment. Adjust before submission!
outfile1 = ['chla_' region '.png'];
print(h, '-dpng', '-r300',  outfile1)
set(gcf, 'renderer', 'painters')
outfile2 = ['chla_' region '.eps'];
print(h, '-depsc2',  outfile2)


%Save information on display
path_display = ['~/data/matlab_bathym/chla_maps/' region '_display.mat'];

save(path_display, 'paperwidth', 'paperheight', 'cbar_width', ...
     'cbar_offset', 'offset2', 'tipos1', 'cbar_loc', 'north', 'lefs', ...
     'rigs', 'tops', 'bots', 'figw', 'figh');
