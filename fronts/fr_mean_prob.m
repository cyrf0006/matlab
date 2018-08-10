function fr_mean_prob(sstList, LatLon, years, months, varargin)

% function fr_mean_prob(sstList, LatLon, years, months, destination, varargin)
%
% Usage ex:
%   fr_mean_prob('SST_pacific', '~/data/front_data/PacificLatLon.mat', [1986:2010], [1:12], 'allMean')
%
%    
% !ls /media/Seagate1TB/Front_data/SST_atlantic/*.mat > SST_atlantic
% !ls /media/Seagate1TB/Front_data/SST_atlantic/*.mat > SST_atlantic
%
% Outfile name will be with varargin extension, in local folder. complete this..
% Shoul be run in ~IML/Fronts/matlab_workspace/probability/ and
% then mv into OUTPUT


% Deal with vararagin
if isempty(varargin)==1 % default
    meanType = 'allMean';
elseif size(varargin,2)==1
    meanType = varargin{1};
else
    disp('Wrong input... try "fr_mean_prob"')
    return
end  


% $$$ % Input that should be passed as arguments:
% $$$ 
% $$$ %sstList = 'SST_atlantic';
% $$$ %sstList = 'SST_hudson';
% $$$ %sstList = 'SST_baffin';
% $$$ sstList = 'SST_pacific';
% $$$ 
% $$$ % $$$ LatLon = '~/data/front_data/AtlanticLatLon.mat';
% $$$ % $$$ LatLon = '~/data/front_data/BaffinLatLon.mat';
% $$$ % $$$ LatLon = '~/data/front_data/HudsonLatLon.mat';
% $$$ LatLon = '~/data/front_data/PacificLatLon.mat';
% $$$ 
dateList = '/tmp/tmp';
% $$$ years = [1986:2010];
% $$$ months = [1:12];


%destination = ['~IML/Fronts/matlab_workspace/probability/OUTPUT'];
outname = ['./' sstList '_' meanType '_prob.mat'];
% Maybe I could keep meanType = [] for allMean...


% I should have a varargin that could be 'yearly', 'monthly', 'season' 
yearly = 1;
monthly = 0;
seasonaly = 0;

fid = fopen(sstList);
C = textscan(fid, '%s', 'delimiter', '\n');
sstFiles = char(C{1});
nFiles = size(sstFiles, 1); %number of eps_files 


% Extrat Date Vector
command = sprintf('sed "s/\\(^.*\\)\\([0-9]\\{8\\}\\)\\(.*$\\)/\\2/" %s > /tmp/tmp', sstList);
system(command);
fid = fopen(dateList);
C = textscan(fid, '%s', 'delimiter', '\n');
dateVec_raw = char(C{1});
dateVec = datenum(str2num(dateVec_raw(:,1:4)), str2num(dateVec_raw(:,5:6)), str2num(dateVec_raw(:,7:8)));

% Latitude longitude
LatLon = load(LatLon);
lat = LatLon.lat;
lon = LatLon.lon;



% Execute function according to varargin input
if strcmp(meanType, 'allMean') == 1 
    I = find(dateVec >= datenum(years(1),1,1) &  dateVec <= datenum(years(end),12,31));
    [edgeCount, pixelCount, probability, cloudProb] = fr_prob_calculation(sstFiles(I,:));
    
    outname = ['./' sstList '_prob.mat'];
    disp(['save ' outname])
    if size(lat,1) ~= size(probability,1)
        disp(['Latitude not corresponding to prob. size. Sure you ' ...
              'call the good one?'])
        keyboard
    else
        save(outname, 'probability', 'pixelCount', 'edgeCount', 'cloudProb', 'lat', 'lon')    
    end
    
elseif strcmp(meanType, 'yearly') == 1 
    for i = 1:length(years)
        I = find(dateVec >= datenum(years(i),1,1) &  dateVec <= datenum(years(i),12,31));
        [edgeCount, pixelCount, probability, cloudProb] = fr_prob_calculation(sstFiles(I,:));
        
        outname = ['./' sstList '_' datestr(datenum(years(i), 1, 1),10) '_prob.mat'];
        disp(['save ' outname])
        if size(lat,1) ~= size(probability,1)
            disp(['Latitude not corresponding to prob. size. Sure you ' ...
              'call the good one?'])
            keyboard
        else
            save(outname, 'probability', 'pixelCount', 'edgeCount', 'cloudProb', 'lat', 'lon') 
        end
    end
    
elseif strcmp(meanType, 'monthly') == 1 
    for i = 1:length(years)
        for j = 1:length(months)
            I = find(dateVec >= datenum(years(i),months(j),1) &  dateVec <= datenum(years(i),months(j),31));
            if isempty(I) == 1
                disp(['Skip ' datestr(datenum(years(i), months(j), 1),28)]);
                continue
            else
                [edgeCount, pixelCount, probability, cloudProb] = fr_prob_calculation(sstFiles(I,:));

                outname = ['./' sstList '_' datestr(datenum(years(i),months(j), 1),10) '-' datestr(datenum(years(i), months(j), 1),7) '_prob.mat'];
                disp(['save ' outname])
                if size(lat,1) ~= size(probability,1)
                    disp(['Latitude not corresponding to prob. size. Sure you ' ...
                          'call the good one?'])
                    keyboard
                else 
                    save(outname, 'probability', 'pixelCount', 'edgeCount', 'cloudProb', 'lat', 'lon') 
                end 
            end
        end
    end
    
elseif strcmp(meanType, 'season') == 1 
    disp('Not implemented yet!')

elseif strcmp(meanType, 'monthlyClim') == 1 
    disp('Not implemented yet!')

elseif strcmp(meanType, 'seasonClim') == 1 
    disp('Not implemented yet!')

else
    disp('Wrong input... try "fr_mean_prob"')
    return
end



disp('FINISHED! [K]')



