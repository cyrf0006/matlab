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
lefs = 0.07; % very left of figure
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
paperwidth = 24;%cm
paperheight = 20;%cm
FS = 12;

% param for colorbar
cbar_width = 0.02;
cbar_offset = 0.01; % colorbar offset from figure
offset2 = 0.1; 


%% Get global bathym    
path = ['~/data/matlab_bathym/baltic.mat']; % 30-sec.
lims = [50 70 0 30];

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
lon_min=11;
lon_max=20;
lat_min=54;
lat_max=57;

I=find(lat<lat_max & lat>lat_min);
J=find(lon<lon_max & lon>lon_min);
latitude=lat(I);
longitude=lon(J);
bathy=z(I,J);


h = figure('visible', 'off');
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 paperwidth paperheight])
m_proj('mercator','long',[min(longitude) max(longitude)],'lat',[min(latitude) max(latitude)]);
%m_grid('box','fancy')
hold on
load ~/PhD/bathym/gebco64
colormap(gebco);
V=[0:10:100];
m_pcolor(longitude,latitude,abs(bathy)); shading flat;
hold on
[HH, HH] = m_contour(longitude,latitude,abs(bathy), [0:20:100], 'color', 'k');
m_gshhs_h('patch',[1 .9333 .6667]); %coastlines (Beige)                               
xlabel('Longitude', 'FontSize', FS, 'fontweight', 'bold')
ylabel('Latitude', 'FontSize', FS, 'fontweight', 'bold')
caxis([min(V) max(V)])
%set(gca, 'fontsize', FS)

myColor1 = [0.0, 0.55, 0.55];
myColor2 = [1.0, 0.0, 0.5];
white = [1 1 1]*.95;

% few annotations
m_line(15.9879, 55.2494,'marker','p','color',white, 'markerfacecolor', myColor1,'markersize',18);
m_line(16.0969, 54.9455,'marker','p','color',white, 'markerfacecolor', myColor2,'markersize',18);
% $$$ m_line(16.0489, 54.9633,'marker','.','color',[1 0 0], 'markerfacecolor', 'r','markersize',14);
% $$$ m_line([15.9945, 16.197],[54.9803, 54.9106],'color','r','linewi',1,'linestyle', '-'); 
m_text(16.05, 55.3,'S1','color',white,'fontSize',16, 'fontweight', 'bold', 'horizontalAlignment', 'left', 'verticalAlignment', 'bottom');
m_text(16.15, 55,'T1B','color',white,'fontSize',16, 'fontweight', 'bold', 'horizontalAlignment', 'left', 'verticalAlignment', 'bottom');

% $$$ m_line(15.9879, 55.2494,'marker','o','color',[1 1 1], 'markerfacecolor', 'm','markersize',6);
% $$$ m_line(15.9924, 54.9795,'marker','o','color',[1 1 1], 'markerfacecolor', 'm','markersize',6);

m_grid('box','fancy')


adjust_space


%% Large Area map (inset)
lon_min=-14;
lon_max=50;
lat_min=32;
lat_max=70;

I=find(lat1m<lat_max & lat1m>lat_min);
J=find(lon1m<lon_max & lon1m>lon_min);
latitude=lat1m(I);
longitude=lon1m(J);
bathy=z1m(I,J);

a2 = axes('position',[0.06 0.485 0.33 0.33]) ; % inset
m_proj('mercator', 'long',[lon_min lon_max],'lat',[lat_min lat_max]);
hold on
%m_pcolor(longitude,latitude,abs(bathy)); shading flat;
%[HH, HH] = m_contour(longitude,latitude,abs(bathy), [0:100:1000], 'color', 'k');
m_gshhs_i('patch',[1 .9333 .6667]); %coastlines (Beige)                               
caxis([min(V) max(V)])
m_grid('box','fancy', 'yticklabels',[], 'xticklabels',[])

% Rectangle
lon_min=11;
lon_max=20;
lat_min=54;
lat_max=57;
[lomi, lami] = m_ll2xy(lon_min,lat_min);
[loma, lama] = m_ll2xy(lon_max,lat_max);

rectangle('Position', [lomi lami loma-lomi lama-lami],'linewidth', 2, 'edgecolor', 'r') 

cbar_offset = 0;%0.04; 
offset2 =0;% 0.02; % offset between heigth of colorbar 
c = colorbar('location','southoutside');
cb_pos = get(c, 'position');
cb_back = cb_pos;
cb_pos(4) = .02;
cb_pos(1) = .6;%cb_pos(1)+.1*cb_pos(3);
cb_pos(3) = .25;
cb_pos(2) = .26;%cb_pos(2)+cbar_offset;
set(c, 'pos', cb_pos);
set(c, 'fontsize', FS, 'fontweight', 'bold')

ti = ylabel(c,'Depth (m)', 'FontSize', FS, 'fontweight', 'bold');
ti_pos = get(ti, 'position');
set(ti, 'Rotation',0.0);     
clim = get(gca, 'clim');
set(ti, 'position', [(clim(2)-clim(1))./2  2 ti_pos(3)]); 
set(gca, 'fontSize', FS)

%save figure
print('-dpng', '-r300', 'BOR10Map.png')
set(gcf, 'renderer', 'painters'); % vectorial figure
print('-depsc', 'BOR10Map.eps')



