%function cfg = sx_readCfg(cfgFile)

%% Other idea, read the whole structure and store in a cell on
%% which I can search for fields.


cfgFile = '~/research/MIO/seaExplorerData/M78-OSCAHR/MIO_SEA007_78/RAW/raw_fulldata/seapayload.cfg';
   
% open file (read-only)
fid=fopen(cfgFile,'rt');

% read header
% $$$ dheader=fgetl(fid);
% $$$ 
% $$$ % Loop on rest of file
% $$$ ddataMat = [];
% $$$ dtimeVec = [];
EOF = 0;
% $$$ cfg = struct();
% $$$ cfg.dpl = {};
% $$$ cfg.sn = {};
% $$$ cfg.calib = {};
s = struct();
device_count = 1;
while ~feof(fid)
    str=fgetl(fid);  
    if str~=-1                    
       
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
            str = fgetl(fid);
            if strfind(str, 'used=yes')
                str = fgetl(fid);
                I = strfind(str, '=');
                device = str(I+1:end);
                command = sprintf('s.device{%d} = ''%s''', device_count, device);
                eval(command);
                command = sprintf('s.slot{%d} = ''%s''', device_count, slot);
                eval(command);
                device_count = device_count + 1;
            else
                continue
            end
            
            % if in sensor configuration
        elseif ~isempty(regexp(str, '\[*]'))
            inDevice = 1;
            while inDevice
                if 
                    
            keyboard

            
        end
        
    else
        continue % empty line            
    end
        

end
