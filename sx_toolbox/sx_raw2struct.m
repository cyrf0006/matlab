function sx_raw2struct(dplFile)
    
% function sx_raw2struct(dplFile)

% Take glider's raw 'science' and 'pilot' files and merge them into 
% a structure containing all variables ans info on the mission. 
% Was originaly sx_egosx.m, but now I am 'DIYing', avoiding ego
% format. 
%
% The only input is the 'deployment file' containing infos on the
% deployment and where to find the raw files.
%
% usage ex. from ~/research/MIO/seaExplorerData/M78-OSCAHR/data_processing: 
% >> sx_raw2struct('M78.dpl')
%
%
% Frederic.Cyr@mio.osupytheas.fr - Jan. 2016
%
%
% ------------------------------------------------ %

% *NOTE: Ask for an example deployment file .dpl (!)

tic
% Read deployment file
disp(' ------- Header ------ ')
fid = fopen(dplFile);
if fid == -1
    disp('No deployment file was found, please check this! [Quit]')
    return
else
    tline = fgetl(fid);
    while ischar(tline)
        eval(tline);
        tline = fgetl(fid);
    end
    fclose(fid);
end
disp(' ---------------------  ')

% load 2 lists in memory
fid = fopen(dataList);
C = textscan(fid, '%s', 'delimiter', '\n');
dataFiles = char(C{1});

fid = fopen(gliderList);
C = textscan(fid, '%s', 'delimiter', '\n');
gliderFiles = char(C{1});

nodFiles = size(dataFiles, 1);
nogFiles = size(gliderFiles, 1);

if size(dataFiles, 1) ~= size(gliderFiles, 1)
    disp('[WARNING] 2 lists not the same size!')
end

% Initialization
datMatrix = [];
logMatrix = [];
mtimeDat = [];
mtimeLog = [];
datFiles = [];
logFiles = [];


%% --- DATA file --- %%
disp('-> Science Files')
for i = 1:nodFiles
    dfile = dataFiles(i, :);
    % Just in case there are white space in filename
    I = find(dfile==' '); dfile(I) = [];

    % open file (read-only)
    fid=fopen(dfile,'rt');

    % read header
    dheader=fgetl(fid);
    
    % Loop on rest of file
    ddataMat = [];
    dtimeVec = [];
    EOF = 0;
    while ~EOF
        str=fgetl(fid);
        if str~=-1            
            I = find(str==';');
            theDate = str(1:I(1)-1);
            
            theVector = [];
            for ii = 1:length(I)-1
                element = str2num(str(I(ii):I(ii+1)));
                if isempty(element)
                    element = NaN;
                end
                theVector = [theVector element];
            end
            
            ddataMat = [ddataMat; theVector];
            dtimeVec = [dtimeVec datenum(theDate, 'dd/mm/yyyy HH:MM:SS.FFF')];
        else
            EOF = 1; % end-of-file
        end
    end    
    
    if i == 1 % variable names
        datList = strread(dheader,'%s','delimiter',';');
        datList(1) = []; % remove 'timestamp'

        % *In early SeaEx deploy, variables names were messed.
        %   **for MFL serial no in [5->8], V3 and V4 inverted!!
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
    I = find(dfile=='/'); dfilename = dfile(I(end)+1:end);
    disp(dfilename);
    datFiles{i} = dfilename;
end
    
%% --- GLIDER file --- %%
disp('-> Navigation Files')
for i = 1:nogFiles
    gfile = gliderFiles(i, :);
    I = find(gfile==' '); gfile(I) = [];    

    % open file (read-only)
    fid=fopen(gfile,'rt');

    % read header
    gheader=fgetl(fid);
    
    % Loop on rest of file
    gdataMat = [];
    gtimeVec = [];
    EOF = 0;
    while ~EOF
        str=fgetl(fid);
        if str~=-1            
            I = find(str==';');
            theDate = str(1:I(1)-1);
            
            theVector = [];
            for ii = 1:length(I)-1
                element = str2num(str(I(ii):I(ii+1)));
                if isempty(element)
                    element = NaN;
                end
                theVector = [theVector element];
            end
                        
            gdataMat = [gdataMat; theVector];
            gtimeVec = [gtimeVec datenum(theDate, 'dd/mm/yyyy HH:MM:SS')];
        else
            EOF = 1; % end-of-file
        end
    end    
    
    if i == 1 % variable names
        logList = strread(gheader,'%s','delimiter',';');
        logList(1) = [];
    end

    % Concatenate nav data
    logMatrix = [logMatrix; gdataMat];
    mtimeLog = [mtimeLog gtimeVec];
    I = find(gfile=='/'); gfilename = gfile(I(end)+1:end);
    disp(gfilename);
    logFiles{i} = gfilename;
end

%% Fill structure
% Header
s = struct();
s.glider_name = glider_name;
s.obs_name = obs_name;
s.campaign_name = campaign_name;
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
