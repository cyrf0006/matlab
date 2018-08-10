function sx_surveyMap(structFile, varargin)
    
    
% function sx_surveyMap(structFile, [lims])
%  
%  where 'lims' an optional argument and whould correspond to:
%  [lat_min lat_max lon_min lon_max]
%
%  usage ex: sx_surveyMap('osca01_struct.mat', [42 45 6 10])
% 
    
%% Load structure and get glider positions
load(structFile)
w = whos;

found = 0;
i = 1;
while found == 0
    if strcmp(w(i).class, 'struct')
        found = 1;
        structName = w(i).name;
    else
        i = i+1;
    end
end
command = sprintf('s = %s ;',structName);
eval(command);

latDeg = floor(s.log.Lat/100);
lonDeg = floor(s.log.Lon/100);
latMin = s.log.Lat - floor(s.log.Lat/100)*100;
lonMin = s.log.Lon - floor(s.log.Lon/100)*100;
latVec = latDeg + latMin/60;
lonVec = lonDeg + lonMin/60;
timeVec = s.mtimeLog;

I = find(diff(latVec) == 0 & diff(lonVec) == 0); % avoid repetition
latVec(I) = [];
lonVec(I) = [];
timeVec(I) = [];
I  = find(lonVec == 0 | latVec == 0); % remove zeros
latVec(I) = [];
lonVec(I) = [];
timeVec(I) = [];

I = find(latVec>=-90 & latVec<=90 & lonVec>=0 & lonVec<=360);
latVec = latVec(I);
lonVec = lonVec(I);
timeVec = timeVec(I);

%% Deal with varargin (map limits)
if isempty(varargin)
    lonLims = [min(lonVec)-2 max(lonVec)+2];
    latLims = [min(latVec)-2 max(latVec)+2];        
elseif length(varargin{1} == 4)
    lims = varargin{1};
    lonLims = lims(1:2);
    latLims = lims(3:4);
else
    disp('Wrong input, check help menu [Quit]')
    return
end



%% plot map
close all
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.0 ; % horiz. space between subplots
dy = 0.0; % vert. space between subplots
lefs = 0.09; % very left of figure
rigs = 0.25; % very right of figure
tops = 0.02; % top of figure
bots = 0.02; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %
FS = 10;
%V = 100:1:300;

figure(1);
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 20 20])
m_proj('mercator','long',lonLims,'lat',latLims);
hold on
m_gshhs_h('patch',[1 .9333 .6667]); %coastlines (Beige)                               
xlabel('Longitude', 'FontSize', FS, 'fontweight', 'bold')
ylabel('Latitude', 'FontSize', FS, 'fontweight', 'bold')

% add track
h1 = m_plot(lonVec, latVec, '.', 'color', 'm');
% add time
for i = [1 round(length(timeVec)/2) length(timeVec)] % plot one out of X
    h2 = m_text(lonVec(i), latVec(i), datestr(timeVec(i), 31), 'fontWeight', 'bold');
    h3 = m_plot(lonVec(i), latVec(i), '.k', 'markerSize', 8);
    h1 = [h1 h2 h3];
end

m_grid('box','fancy')
adjust_space

% $$$ 
% $$$ set(gcf, 'renderer', 'painters'); % vectorial figure
print('-dpng', '-r300',  [structName '_track.png'])


