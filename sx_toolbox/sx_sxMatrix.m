function sx_sxMatrix(dataList, gliderList, prefix)
    
% usage ex:
% sx_sxMatrix('data.list', 'glider.list', mio001)
    
% load 2 lists in memory
fid = fopen(dataList);
C = textscan(fid, '%s', 'delimiter', '\n');
dataFiles = char(C{1});

fid = fopen(gliderList);
C = textscan(fid, '%s', 'delimiter', '\n');
gliderFiles = char(C{1});

if size(dataFiles, 1) == size(gliderFiles, 1)
    noFiles = size(dataFiles, 1);
else
    disp('[ERROR] 2 lists not the same size!')
    return
end


for i = 1:noFiles
    dfile = dataFiles(i, :);
    gfile = gliderFiles(i, :);
    disp(dfile)
    % Just in case there are white space in filename
    I = find(dfile==' '); dfile(I) = [];
    I = find(gfile==' '); gfile(I) = [];    
    
    
    %% --- DATA file --- %%
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
            dtimeVec = [dtimeVec datenum(theDate, 'dd/mm/yyyy HH:MM:SS')];
        else
            EOF = 1; % end-of-file
        end
    end
    
        
    %% --- GLIDER file --- %%
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
    
    %% Save files
    outFileScience = [prefix '_science_' sprintf('%.04d', i) '.mat'];
    outFilePilot = [prefix '_pilot_' sprintf('%.04d', i) '.mat'];
    
    dvariables = strread(dheader,'%s','delimiter',';') 
    gvariables = strread(gheader,'%s','delimiter',';') 
    dvariables(1) = []; % remove 'timestamp'
    gvariables(1) = [];

    % science
    data = [];
    data.matrix = ddataMat;
    data.mtime = dtimeVec;
    data.name = dvariables;
    save(outFileScience, 'data')
    
    % pilot
    data = [];
    data.matrix = gdataMat;
    data.mtime = gtimeVec;
    data.name = gvariables;
    save(outFilePilot, 'data')    
    
    keyboard
end

