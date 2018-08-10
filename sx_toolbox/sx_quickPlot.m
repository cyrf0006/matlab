function sx_quickplot(structFile, varargin)
    
% function sx_varnames(gliderFile)
%
% Function use to quickly check variables names in a structure file 
% 
% usage ex:
%  sx_quickplot('osca01_struct.mat', 'SBD_TEMPERATURE', 'SBD_CONDUCTIVITY')
%  sx_quickplot('osca01_struct.mat', 'UV1_PHE', 'UV1_TRY')
    
%deal with varargin    
if isempty(varargin)
    disp(['!!! you should ask for at least one variable !!! (see sx_varnames.m ' ...
          'for names)'])
    return
else
    
    % get data
    load(structFile);
    w = whos;
    structName = w.name;
    command = sprintf('s = %s ;',structName);
    eval(command);    
    
    z = s.data.NAV_DEPTH;
    n = s.mtimeDat;       
        
    % interp z to time (remove NaNs)
    I = find(~isnan(z));
    z = interp1(n(I), z(I), n);
        
    % "dirty" cast detection
    W = gradient(z)./gradient(n*86400); % m/s
    dt = diff(abs(n(1:2)))*86400; %sec
    fs = 1/dt; % Hz
    freq_low = 0.01; %Hz
    Wn_low = freq_low/(fs/2);
    [b,a] = butter(4, Wn_low);
    I= find(isnan(W));
    W(I) = 0;
    Wfilt = filtfilt(b, a, W);

    I = find(abs(Wfilt)<.01);
    z(I) = NaN; % flag as NaN when not moving
    
    % build the grid
    x = min(n):5/1440:n(end);
    y = 0:max(z);
    [X,Z] = meshgrid(x,y);    
    
    for i = 1:length(varargin)
        disp(sprintf('Evaluating variable: %s', varargin{i})); 
        
        command = sprintf('myVar = s.data.%s ;', varargin{i});
        eval(command);
                
        % grid data
        I = find(~isnan(myVar) & ~isnan(z'));
        [XI,ZI,VAR] = griddata(n(I),z(I),myVar(I),X,Z);
        
        figure(i)
        clf
        imagesc(x,y, VAR)
        datetick
        colorbar
        title(varargin{i})
        ylabel('Depth (m)')
        xlim([min(x) max(x)])
        
        if i<length(varargin)
        disp(' Press any key for next plot')
        pause
        end
        
    end
end

disp('Done!')
