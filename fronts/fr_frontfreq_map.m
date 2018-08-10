function fr_frontfreq_map(frontData, region, varargin)

%usage ex:    
% fr_frontfreq_map('~/research/fronts/matlab_workspace/probability/OUTPUT/SST_atlantic_prob.mat','Atlantic', [30 66 -77 -35])
% fr_frontfreq_map('~/research/fronts/matlab_workspace/probability/OUTPUT/SST_atlantic_prob.mat', 'NSshelf', [41 45 -72 -58])
% fr_frontfreq_map('~/research/fronts/matlab_workspace/probability/OUTPUT/SST_atlantic_prob.mat', 'LabradorShelf', [52 65 -79 -54], 'East')
% fr_frontfreq_map('~/research/fronts/matlab_workspace/probability/OUTPUT/SST_atlantic_prob.mat', 'NFDLinner', [44 52 -60 -51], 'East')    
% fr_frontfreq_map('~/research/fronts/matlab_workspace/probability/OUTPUT/SST_atlantic_prob.mat', 'NFDLshelf', [42 54 -60 -43], 'East')    
% fr_frontfreq_map('~/research/fronts/matlab_workspace/probability/OUTPUT/SST_atlantic_prob.mat', 'NFDLshelf', [40 47 -60 -44])
% fr_frontfreq_map('~/research/fronts/matlab_workspace/probability/OUTPUT/SST_atlantic_prob.mat', 'HudsonSt', [58 65 -79 -60], 'North')
% fr_frontfreq_map('~/research/fronts/matlab_workspace/probability/OUTPUT/SST_atlantic_prob.mat', 'SouthLabrador', [50 59 -63 -54], 'East')
% fr_frontfreq_map('~/research/fronts/matlab_workspace/probability/OUTPUT/SST_atlantic_prob.mat', 'GSL', [45 52 -70 -55])

% fr_frontfreq_map('~/research/fronts/matlab_workspace/probability/OUTPUT/SST_hudson_prob.mat', 'Hudson', [50 71 -96 -72], 'East')
% fr_frontfreq_map('~/research/fronts/matlab_workspace/probability/OUTPUT/SST_hudson_prob.mat', 'James', [51 56 -83 -78], 'East')
% fr_frontfreq_map('~/research/fronts/matlab_workspace/probability/OUTPUT/SST_hudson_prob.mat', 'HudsonSouth', [55 60 -96 -84], 'North')
% fr_frontfreq_map('~/research/fronts/matlab_workspace/probability/OUTPUT/SST_hudson_prob.mat', 'HudsonWest', [60 64 -95 -90], 'East')
% fr_frontfreq_map('~/research/fronts/matlab_workspace/probability/OUTPUT/SST_hudson_prob.mat', 'Belcher', [55 59 -84 -76], 'East')
% fr_frontfreq_map('~/research/fronts/matlab_workspace/probability/OUTPUT/SST_hudson_prob.mat', 'Foxe', [64 71 -85 -70], 'East')

% fr_frontfreq_map('~/research/fronts/matlab_workspace/probability/OUTPUT/SST_baffin_prob.mat', 'Baffin', [60 82 -85 -50], 'East')

% fr_frontfreq_map('~/research/fronts/matlab_workspace/probability/OUTPUT/SST_pacific_prob.mat', 'Pacific', [30 67 -179 -115], 'North')
% fr_frontfreq_map('~/research/fronts/matlab_workspace/probability/OUTPUT/SST_pacific_prob.mat', 'BeringSea', [50 67 -179 -156], 'East')
% fr_frontfreq_map('~/research/fronts/matlab_workspace/probability/OUTPUT/SST_pacific_prob.mat', 'GulfAlaska', [50 62 -157 -126], 'North')
% fr_frontfreq_map('~/research/fronts/matlab_workspace/probability/OUTPUT/SST_pacific_prob.mat', 'WestCoast', [30 52 -130 -115], 'East')
% fr_frontfreq_map('~/research/fronts/matlab_workspace/probability/OUTPUT/SST_pacific_prob.mat', 'Vancouver', [44 55 -132 -122])


% !!!! IMPORTANT: if you want to used saved display value for a certain
% region, you must only provide 1 varargin argument (otherwise
% display info will be overwritten) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.05; % very right of figure
tops = 0.05; % top of figure
bots = 0.05; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %


%     Preamble, check input and deal with varargin    
path = ['~/data/matlab_bathym/' region '.mat'];

fid = fopen(path);
if fid == -1 % doesnt exist yet!
    disp('Region doesn''t exist yet!')
    if isempty(varargin) == 1 
        disp([' -> [ERROR!] No coordinates provided. Please check help menu'])
        return
    else
        if length(varargin{1}) == 4
            disp(' -> Extract bathym from Gebco 30'''', may take some time...')
            lims = varargin{1};       
            [lat lon z] = getGebco('~/data/GEBCO/GEBCO_08.nc', 'z', lims);
            disp('   done!')
            save(path, 'lat', 'lon', 'z');
        else
            disp(' -> [ERROR!] Wrong input, please check help menu')
            return
        end
    end
else
    load(path)
end



% Check for figure dimensions, colorbar position, etc.
if size(varargin,2) == 2 
    if strcmp(varargin{2}, 'North') == 1
        paperwidth = 23;%cm
        paperheight = 18;%cm
        cbar_width = 0.02;
        cbar_offset = 0.1; % colorbar offset from figure
        offset2 = 0.1; % offset between heigth of colorbar 
        tipos1 = 0.15;
        lefs = 0.1; % very left of figure
        rigs = 0.05; % very right of figure
        tops = 0.05; % top of figure
        bots = 0.05; % bottom of figure
        figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
        figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
        cbar_loc = 'NorthOutside';
        north = 1;
    else        
        % (East)
        paperwidth = 15;%cm
        paperheight = 25;%cm
        cbar_width = 0.02;
        cbar_offset = 0.12; % colorbar offset from figure
        offset2 = 0.1; % offset between heigth of colorbar 
        tipos1 = 0.15;
        cbar_loc = 'EastOutside';
        lefs = 0.05; % very left of figure
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
    path_display = ['~/data/matlab_bathym/' region '_display.mat'];
    load(path_display)
end

path_display = ['~/data/matlab_bathym/' region '_display.mat'];
load(path_display)


load('~/data/matlab_bathym/GebcoColormapFred.mat')

% data reduction
front = load(frontData);
frontlat = front.lat(:,1);
frontlon = front.lon(1,:);
prob = front.probability;

I = find(frontlat > lat(end) & frontlat < lat(1));
J = find(frontlon > lon(1) & frontlon < lon(end));
frontlat = frontlat(I);
frontlon = frontlon(J);
prob = prob(I,J);

% $$$ % Smooth bathymetry
% $$$ b=gaussfir(.01);  %Gaussian filter designer
% $$$ z = filter2(b,z);  %Application of the filter
% Probability with z = -1000m isobat
%keyboard
% figure
h = figure('Visible', 'off');
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 paperwidth paperheight])
%set(gcf, 'renderer', 'opengl')
load /home/cyrf0006/research/fronts/matlab_workspace/GS-SB_position.mat
m_proj('mercator','long',[min(lon) max(lon)],'lat',[min(lat) max(lat)]);
m_grid('box','fancy')
hold on
%  WATCH OUT! make sure m_gshhs_h is  'high res' for final figure
m_pcolor(frontlon,frontlat,prob); shading flat;
%m_gshhs_h('patch',[.7 .7 .7]); %coastlines (gray)                           
m_gshhs_h('patch',[1 .9333 .6667]); %coastlines (Beige)                               
%[HH, HH] = m_contour(lon,lat,z, [0 0], 'color', 'k');
%[HH, HH] = m_contour(lon,lat,z, [0 -100 -200 -300 -500 -1000], 'color', 'k');
%[H2, H2] = m_contour(lon,lat,z, [-1000 -1000], 'color', 'k');
%[H3, H3] = m_contour(lon,lat,z, [-25 -50 -100], 'color', 'k');
% $$$ m_plot(ShelfBreak(:,1), ShelfBreak(:,2), '--r', 'linewidth', 2)
% $$$ m_plot(GulfStream(:,1), GulfStream(:,2), 'r', 'linewidth', 2) 

plot_border('~/research/fronts/borders/NorthAmerica.dat',[min(lon) max(lon)], [min(lat) max(lat)]);



% ---------- Xtra features on different plots ------------- %
%      (see also external function fr_putplaces.m)
if strcmp(region, 'HudsonSt') == 1
    load ~/research/fronts/matlab_workspace/taggart_line.txt;
    taggart_line(end,:) = [];
    m_plot([taggart_line(:,1)], [taggart_line(:,2)], 'r', 'linewidth', 2);
    [HH, HH] = m_contour(lon,lat,z, [-50 -200 -1000], 'color', 'k');
    fr_putPlaces('HUDS')
end

if strcmp(region, 'SouthLabrador') == 1
    [HH, HH] = m_contour(lon,lat,z, [-50 -200 -1000], 'color', 'k');
    fr_putPlaces('LABS')
end

if  strcmp(region, 'NSshelf') == 1
    [HH, HH] = m_contour(lon,lat,z, [-50 -100 -200 -1000], 'color', 'k');
    fr_putPlaces('SCOS')
end

if strcmp(region, 'NFDLinner') == 1 | strcmp(region, 'NFDLshelf') == 1
    [HH, HH] = m_contour(lon,lat,z, [-100 -200 -1000], 'color', 'k');
    fr_putPlaces('NFDL')
end

if strcmp(region, 'Baffin') == 1
    [H3, H3] = m_contour(lon,lat,z, [-500 -1000], 'color', 'k');
    load ~/research/fronts/matlab_workspace/lobb2003.mat
    fr_putPlaces('BAFF')
    %m_plot(lon(3:4), lat(3:4), 'm', 'linewidth', 2)
    m_line(-69.3457, 77.4830, 'marker','.','MarkerFaceColor', ...
           'k','markersize',14,'color',[.6 0 0]);
end

if strcmp(region, 'GSL') == 1
    [HH, HH] = m_contour(lon,lat,z, [-50 -300], 'color', 'k');
    fr_putPlaces('GULF')
    m_line(-69.716606, 48.156604,'marker','.','MarkerFaceColor', ...
           'k','markersize',14,'color',[.6 0 0]);
    m_line(-67.386367, 49.35, 'marker','.','MarkerFaceColor', ...
           'k','markersize',14,'color',[.6 0 0]);
end

if strcmp(region, 'Pacific') == 1
    [H3, H3] = m_contour(lon,lat,z, [-1000 -1000], 'color', 'k');
    fr_putPlaces('PACI')
end

if strcmp(region, 'Atlantic') == 1
    [H3, H3] = m_contour(lon,lat,z, [-1000 -1000], 'color', 'k');
end

if strcmp(region, 'Hudson') == 1 
    [H3, H3] = m_contour(lon,lat,z, [-100 -100], 'color', 'k');
    HudsonRiver = load('~/research/fronts/rivers/rivers_clim_lat_lon');
    %load ~/LaTeX/MS/fronts/HudsonRiver.dat;
    %    indices = [4 5 7 8 10:23];
    indices = [1:7 27 28 37:39 45:50 60];
% $$$     for i = 1:length(indices)
% $$$         lat = HudsonRiver(indices(i), 3);
% $$$         lon = HudsonRiver(indices(i), 4);
% $$$         m_line(lon, lat,'marker','o','MarkerFaceColor','m', ...
% $$$                'markersize',5,'color',[1 1 1]);
% $$$     end
    fr_putPlaces('HUDB')
end

if   strcmp(region, 'Belcher') == 1
    [H3, H3] = m_contour(lon,lat,z, [-25 -50 -100], 'color', 'k');
end

if   strcmp(region, 'HudsonWest') == 1
    [H3, H3] = m_contour(lon,lat,z, [-25 -50 -100], 'color', 'k');
    fr_putPlaces('WHUD')
end

if   strcmp(region, 'Foxe') == 1
    [H3, H3] = m_contour(lon,lat,z, [-25 -50 -100], 'color', 'k');
    fr_putPlaces('FOXE')
end

if   strcmp(region, 'Vancouver') == 1
    [H3, H3] = m_contour(lon,lat,z, [-100 -1000], 'color', 'k');
    fr_putPlaces('VANC')
end

% $$$ if strcmp(region, 'James') == 1 | strcmp(region, 'HudsonSouth') == ...
% $$$         1 | strcmp(region, 'HudsonWest') == 1 | strcmp(region, 'Belcher') == 1 | strcmp(region, 'Foxe') == 1
% $$$     [H3, H3] = m_contour(lon,lat,z, [-25 -50 -100], 'color', 'k');
% $$$     HudsonRiver = load('~/research/fronts/rivers/rivers_clim_lat_lon');
% $$$     %    load ~/LaTeX/MS/fronts/HudsonRiver.dat;
% $$$     %indices = [4 5 7 8 10:23];
% $$$     indices = [1:7 27 28 37:39 45:50 60];
% $$$ % $$$     for i = 1:length(indices)
% $$$ % $$$         lat = HudsonRiver(indices(i), 3);
% $$$ % $$$         lon = HudsonRiver(indices(i), 4);
% $$$ % $$$         % m_text
% $$$ % $$$         m_line(lon, lat,'marker','o','MarkerFaceColor','m', ...
% $$$ % $$$                'markersize',5,'color',[1 1 1]);
% $$$ % $$$     end
% $$$ end

if strcmp(region, 'BeringSea') == 1
    [H3, H3] = m_contour(lon,lat,z, [-50 -100 -1000], 'color', 'k');
    m_line([-172 -167.5], [54 55], 'color', 'k','linewidth', 2);
    m_line([-173 -166], [56 56.5], 'color', 'k','linewidth', 2);
    m_line([-175 -167], [57.5 58.5], 'color', 'k','linewidth', 2);
    m_text(-172, 54, 'Shelf Break ', 'horizontalAlignment', 'right', 'fontWeight', 'bold')
    m_text(-173, 56, 'Middle ', 'horizontalAlignment', 'right', 'fontWeight', 'bold')
    m_text(-175, 57.5, 'Inner ', 'horizontalAlignment', 'right', 'fontWeight', 'bold')
end

if strcmp(region, 'GulfAlaska') == 1 
    [H3, H3] = m_contour(lon,lat,z, [-100 -1000], 'color', 'k');
    fr_putPlaces('GALA')
end

if  strcmp(region, 'WestCoast') == 1
    [H3, H3] = m_contour(lon,lat,z, [-100 -1000], 'color', 'k');
end
% --------------------------------------------------------- %
hold off

xlabel('Longitude', 'FontSize', 10)
ylabel('Latitude', 'FontSize', 10)
set(gca, 'fontsize', 10)

adjust_space

caxis([0 12])
%caxis([.5 1.5])

m_grid('box','fancy')

if north
    c = colorbar(cbar_loc);
    Pos = get(gca, 'position');
    set(c, 'FontSize', 10, 'position', [Pos(1)+offset2 Pos(2)+Pos(4)+cbar_offset Pos(3)-2*offset2 cbar_width]);
    ti = ylabel(c,'f(%)', 'FontSize', 10);
    ti_pos = get(ti, 'position');
    set(ti, 'Rotation',0.0);     
    clim = get(gca, 'clim');
    set(ti, 'position', [(clim(2)-clim(1))./2  8 ti_pos(3)]); 

else
    c = colorbar(cbar_loc);
    Pos = get(gca, 'position');
    set(c, 'FontSize', 10, 'position', [Pos(1)+Pos(3)+cbar_offset Pos(2)+offset2 cbar_width Pos(4)-2*offset2]);
    ti = ylabel(c,'f(%)', 'FontSize', 10);
    ti_pos = get(ti, 'position');
    set(ti, 'position', [ti_pos(1)-tipos1 ti_pos(2) ti_pos(3)]); 
end
%  WATCH OUT! Here -r100 is quite low resolution, just to speed up
%  treatment. Adjust before submission!
set(gcf, 'renderer', 'painters')
outfile1 = ['freq_' region '.png'];
print(h, '-dpng', '-r300',  outfile1)
% $$$ outfile2 = ['freq_' region '.eps'];
% $$$ print(h, '-depsc2',  outfile2)



% Save information on display
path_display = ['~/data/matlab_bathym/' region '_display.mat'];

save(path_display, 'paperwidth', 'paperheight', 'cbar_width', ...
     'cbar_offset', 'offset2', 'tipos1', 'cbar_loc', 'north', 'lefs', ...
     'rigs', 'tops', 'bots', 'figw', 'figh');




% Cold anomaly position (for Fred's 3rd chapter):
% $$$ lon_min=-69.5452;
% $$$ lon_max=-69.4702;
% $$$ lat_min=48.1076;
% $$$ lat_max=48.1576;
% $$$ [lomi, lami] = m_ll2xy(lon_min,lat_min);
% $$$ [loma, lama] = m_ll2xy(lon_max,lat_max);
% $$$ 
% $$$ rectangle('Position', [lomi lami loma-lomi lama-lami], 'linewidth', ...
% $$$           2, 'edgecolor', 'r') 
