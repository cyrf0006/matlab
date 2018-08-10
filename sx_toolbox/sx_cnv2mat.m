function s = sx_cnv2mat(cnvList)
    
% function sx_cvn2mat(cnvList)

% (TO BE EDITED) Take glider's raw 'science' and 'pilot' files and merge them into 
% a structure containing all variables ans info on the mission. 
% Was originaly sx_egosx.m, but now I am 'DIYing', avoiding ego
% format. 
%
% The only input is the 'deployment file' containing infos on the
% deployment and where to find the raw files.
%
% usage ex. from ~/research/MIO/seaExplorerData/MM230-ROMARIN/data_processing:
% >> sx_cvn2mat('MM230cnv.list')
%
%
% Frederic.Cyr@mio.osupytheas.fr - Feb. 2016
%
%
% ------------------------------------------------ %

% load list in memory
fid = fopen(cnvList);
C = textscan(fid, '%s', 'delimiter', '\n');
cnvFiles = char(C{1});

% Empty structure to be filled
s = struct();

%% Read 1st file to get variable names
file = cnvFiles(1, :);
fid=fopen(file,'rt'); % read-only

inHeader = 1;
varCount = 1;
while inHeader 
    tline  = fgetl(fid);

    if isfloat(tline(1))
        inHeader = 0;
    else
        if strcmp(tline(1), '#') %variable
            I = findstr(tline, ' ');
            tline(I) = []; % remove white spaces
            I = findstr(tline, '='); % position of '='
            J = findstr(tline, ':');
            K = findstr(tline, '[');
            L = findstr(tline, ']');
            
            vshort = tline(I+1:J-1);
            vlong = tline(J+1:K-1);
            units = tline(K+1:L-1);
            
            s.variables.shortNames(varCount) = {vshort};
            s.variables.longNames(varCount) = {vlong};
            s.variables.unitsNames(varCount) = {units};
            varCount = varCount + 1;
        end
    end
end


%% Loop on all files
for i = 1:size(cnvFiles,1)
    file = cnvFiles(i, :);
 
    % Just in case there are white space in filename
    I = find(file==' '); file(I) = [];
    
    % read file line-by-line and fill struct
    fid=fopen(file,'rt'); % read-only

    % Read header
    inHeader = 1; 
    while inHeader
        tline  = fgetl(fid);
        if strcmp(tline(1), '*') %variable
            I = findstr(tline, '='); % position of '='
            if ~isempty(findstr(tline, 'Time'));
                n = datenum(tline(I+2:end), 'mmm dd yyyy HH:MM:SS');
            elseif ~isempty(findstr(tline, 'Latitude'));
                lat = str2num(tline(I+2:end-2));
            elseif ~isempty(findstr(tline, 'Longitude'));
                lon = str2num(tline(I+2:end-2));
            else
                disp('WARNING: some info not read from CVN file')
            end
        elseif strcmp(tline(1), '#') %variable
            continue
        else % Now in data
            inHeader = 0;
            dataMatrix = str2num(tline); % 1st matrix line
        end
    end
    
    % fill dimensions
    s.dimensions.mtime(i) = n;
    s.dimensions.latitude(i) = lat(1) + lat(2)/60;
    s.dimensions.longitude(i) = lon(1) + lon(2)/60;
    
    % Read data
    EOF = 0;
    while ~EOF
        tline  = fgetl(fid);
        if tline == -1
            EOF = 1;
        else
            dataMatrix = [dataMatrix; str2num(tline)];
        end
    end
    
    
    if size(dataMatrix, 2) ~= length(s.variables.shortNames)
        disp(['!!!Problem in files, number of columns do not match number of variables!!!'])
        return
    end
        
    %%%%%%%% HERE!!!!!!!!!!!!!!!!! %%%%%%%%%%%%%
    
    if i == 1 % variable names
        datList = strread(dheader,'%s','delimiter',';');
        logList = strread(gheader,'%s','delimiter',';');
        datList(1) = []; % remove 'timestamp'
        logList(1) = [];

        % *In early SeaEx deploy, variables names were messed.
        %   **for MFL serial no in [5->8], V3 and V4 inverted!!
        keyboard
        for ii = 1:length(datList)
            varName = datList(ii);
            if strcmp(varName, 'UV1_TMP1');
                datList{ii} = 'MFL_TMPE';
            elseif strcmp(varName, 'UV1_TMP2');
                datList{ii} = 'MFL_TMPD';
            elseif strcmp(varName, 'UV1_PHE') | strcmp(varName, 'UV1_V1');
                datList{ii} = 'MFL_V1';
            elseif strcmp(varName, 'UV1_TRY') | strcmp(varName, 'UV1_V2');
                datList{ii} = 'MFL_V2';
            elseif strcmp(varName, 'UV1_LD1') | strcmp(varName, 'UV1_V3');
                datList{ii} = 'MFL_V3';                
            elseif strcmp(varName, 'UV1_LD2') | strcmp(varName, 'UV1_V4');
                datList{ii} = 'MFL_V4';  
            end
        end
        
    end
    
    % Concatenate science data
    datMatrix = [datMatrix; ddataMat];
    mtimeDat = [mtimeDat dtimeVec];
        
    % Concatenate pilot data
    logMatrix = [logMatrix; gdataMat];
    mtimeLog = [mtimeLog gtimeVec];
    
    % Other info
    I = find(dfile=='/'); dfilename = dfile(I(end)+1:end);
    I = find(gfile=='/'); gfilename = gfile(I(end)+1:end);
    disp(dfilename);
    datFiles{i} = dfilename;
    logFiles{i} = gfilename;
end


%% Fill structure
% Header
s = struct;
s.glider_name = glider_name;
s.glider_id = glider_id;
s.obs_name = obs_name;
s.chief_scientist = chief_scientist;
s.glider_respo = glider_respo;
s.deployment_file = dplFile;
s.creation_date = date;

% sensors info
w = [whos('SN*'); whos('calib*')];
for i = 1:length(w)
    sensorName = w(i).name;
    command = sprintf('s.sensors.%s = eval(sensorName) ;', sensorName);
    eval(command);    
end


% dimensions
s.mtimeDat = mtimeDat;
s.mtimeLog = mtimeLog;
s.datVarNames = datList;
s.logVarNames = logList;

% Science data
for i = 1:length(datList)
    command = sprintf('s.data.%s = datMatrix(:,i) ;', datList{i});
    eval(command);
end

% Log/pilot data
for i = 1:length(logList)
    command = sprintf('s.log.%s = logMatrix(:,i) ;', logList{i});
    eval(command);
end

% Info on files:
s.datFiles = datFiles;
s.logFiles = logFiles;


% just rename structure
command = sprintf('%s = s', obs_name);
eval(command);

%% Save files
outFile = [obs_name '_struct.mat'];
save(outFile, obs_name)
toc
disp(' -> FINISHED!')
