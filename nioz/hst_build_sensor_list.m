function theList = hst_build_sensor_list(sensor_id, fileMinSize, varargin)
    
% function sensorList = hst_build_sensor_list(sensor_id, varargin)
%
% This function create a list (with full path) for the location of
% thermistors file considering their serial numbers and folder(s)
% where they can be found.
%
% usage ex (from Rockall Bank 2012):  
% >> hst_build_sensor_list(sensor_id, '/media/Seagate1TB/NIOZ/thermistordata/ROC12/mod4c', '/media/Seagate1TB/NIOZ/thermistordata/ROC12/oldmod4') 
%    OR
% >> hst_build_sensor_list(sensor_id, './mod4c', './oldmod4')

    
% F. Cyr - February 2014
        
% remember original path
%mypath = pwd;

% Empty list
%sensorList = struct();

% possible extensions (edit here if new ones are created)
exts = ['.bin'; '.f7f'];

files = [];
%fullpath = [];
for i = 1:size(varargin,2)    
    thePath = varargin{i}; % check if '/' is there
    if thePath(end)~='/'
        thePath = [thePath '/'];
    end
    
    for iext = 1:size(exts,1) % get all files in folders
        dirpath = [thePath '*' exts(iext,:)];
        tmp = dir(dirpath);  

        if ~isempty(tmp) == 1
            for i = 1:length(tmp)
                tmp(i).path = thePath;
            end
            files = [files; tmp];
        end
    end
end

% We have 3 lists: sensor, path, size
sensorList = char(files.name);
pathList = char(files.path);
for i = 1:length(files)
    sizeList(i) = files(i).bytes;
end

% Keep only sensors corresponding to SN.
theList = []; % final list
missCount = 0;

for i = 1:length(sensor_id)
    sn = sensor_id(i);
    I = find(str2num(sensorList(:,1:4)) == sn); % this is weak.
    if ~isempty(I) == 1       

        [Y, J] = max(sizeList(I)); % keep maximum size file  

        if Y < fileMinSize
            theList = [theList; {'missing_too_short'}]; % for missing sensors
            missCount = missCount+1;
            continue
        end
        
        % remove white spaces
        tmpPath = pathList(I(J),:);
        II = find(tmpPath==' ');
        tmpPath(II) = [];
        
        tmpSensor = sensorList(I(J),:);
        II = find(tmpSensor==' ');
        tmpSensor(II) = [];
        
        % Concatenate full path
        theList = [theList; {[tmpPath tmpSensor]}];
    else 
        theList = [theList; {'missing'}]; % for missing sensors
        missCount = missCount+1;
    end
   
end
% cell2char:
%theList = char(theList);
missCount;


filename = 'sensor.list';
fid = fopen(filename, 'w');
for i=1:size(theList,1)
    fprintf(fid, '%s\n', theList{i});
end
fclose(fid);

% cell2char:
%theList = char(theList);
