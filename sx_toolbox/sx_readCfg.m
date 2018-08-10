function header = sx_readCfg(cfgFile)

% function header = sx_readCfg(cfgFile)
%
% Function that read a SeaExplorer configuration file and returns a
% structure 'header' that contains information on the deployment,
% campaign, calibration constant, etc. The only input is the
% standard .cfg file that was used for glider's deployment, which
% would have been manually edited (preferably) after/before the
% deployement to add complementary info that is not standardized by Alseamar
%
%
% usage ex: header = sx_readCfg('./seapayload.cfg');
%
% returns:
% header = 
%            dpl_info: [1x1 struct]
%    sensors_serialNo: [1x1 struct]
%       sensors_calib: [1x1 struct]
%       glider_config: [1x1 struct]
%
% Frederic.Cyr@mio.osupytheas.fr - February 2016
% ---------------------------------------------------------- %

%% Dump the whole file in a cell array
fid=fopen(cfgFile,'rt');
configFile = {};
while ~feof(fid)
    str=fgetl(fid); 
    configFile = [configFile; {str}];
end
s.configFile = configFile;


%% Get info on deployment and slots configuration
device_count = 1;
for i = 1:numel(configFile)
    str=configFile{i};
    if ~isempty(str)                    

        % if a comment line
        if str(1) == '#' 
            I = [strfind(str, 'dpl_');... 
                 strfind(str, 'sn_');... 
                 strfind(str, 'calib_');];
            
            if ~isempty(I)
                eval(str(I:end));
            else
                continue
            end
            
            % if global parameter ID
        elseif ~isempty(strfind(str, 'id=')); 
            I = strfind(str, '=');
            id=str(I+1:end);

            % if global parameter datalogger info
        elseif ~isempty(strfind(str, 'datalogger'));
            eval(str);

            % if in slot configuration
        elseif ~isempty(regexp(str, '\[J\d\]'))
            slot = str(2:3);
            i = i+1;
            str=configFile{i};
   
            if ~isempty(strfind(str, 'used=yes'))
                i = i+1;
                str=configFile{i};
                I = strfind(str, '=');
                device = str(I+1:end);
                command = sprintf('s.devices{%d} = ''%s''', device_count, device);
                eval(command);
                command = sprintf('s.slots{%d} = ''%s''', device_count, slot);
                eval(command);
                device_count = device_count + 1;
            else
                continue
            end            
        end
    else
        continue % empty line            
    end
end

%% Get info on devices configuration
for i = 1:numel(s.devices)
    Index = strfind(configFile, sprintf('[%s]', s.devices{i}));
    I = find(not(cellfun('isempty', Index)));
    inConfig = 1;

    while inConfig 
        I = I+1;
        
        if I >= numel(configFile)
            inConfig = 0;
            continue
        end
        
        str=configFile{I};
        if isempty(str)
            continue
        elseif str(1) == '#' % at some point we check for SN reading here
            continue
        elseif str(1) == '['
            inConfig = 0;
        elseif ~isempty([strfind(str, 'cfg'); strfind(str, 'acq')])
            command = [sprintf('s.slotsConfig{%d}.', i) str]
            II = regexp(command, '\.\d\.');
            if ~isempty(II); % funky naming convention .1. (to be changed)
                command = [command(1:II) 'p' command(II+1:end)];
            end
            eval(command)
        end
    end
end


%% Fill header with info
header = struct();

% deployment info (only recognized fields)
if exist('dpl_obs_name')
    header.dpl_info.obs_name = dpl_obs_name;
end
if exist('dpl_campaign_name')
    header.dpl_info.campaign_name = campaign_name;
end
if exist('dpl_chief_scientist')
    header.dpl_info.chief_scientist = dpl_chief_scientist;
end
if exist('dpl_glider_respo')
    header.dpl_info.glider_respo = dpl_glider_respo;
end

% sensors serial no. ("sn_" pattern);
w = whos('sn_*');
for i = 1:numel(w)
    varName = w(i).name;
    command = sprintf('header.sensors_serialNo.%s = %d ;', varName, eval(varName));
    eval(command);
end
% !!!!!!! Might not be needed because in configFile anyway !!!!!
%          (but useful to correct deploy #MLF<8)


% calibration constants ("calib_" pattern)
w = whos('calib_*');
for i = 1:numel(w)
    varName = w(i).name;
    varName
    if isempty(eval(varName))
        command = sprintf('header.sensors_calib.%s = %d ;', varName, NaN);
    elseif ischar(eval(varName))
        command = sprintf('header.sensors_calib.%s = ''%s'' ;', varName, eval(varName));
    elseif isnumeric(eval(varName))
        command = sprintf('header.sensors_calib.%s = %d ;', varName, eval(varName));
    end
    eval(command);    
end

% Put structure 's' (glider's config) in header
header.glider_config = s;
