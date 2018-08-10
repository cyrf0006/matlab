function hst_h52mat(hdfFile, nSkip, matFile, varargin)

% usage ex: 
%      hst_h52mat('BOR1001_S1_p4std.h5', 0, 'Titp_S1_p4std.mat', 'S1')
%      hst_h52mat('BOR1001_T1B_p2std.h5', 0, 'Titp_T1B_p2std.mat', 'T1B')
%      hst_h52mat('BOR1001_T1C_p4std.h5', 0, 'Titp_T1C_p4std.mat', 'T1C')
%      hst_h52mat('myFilename.h5', 0, 'rockall_full.mat')
%      hst_h52mat('BOR1001_T1C_reco_4std.h5', 0, 'Titp_T1C_p4std.mat', 'T1C')



% Must adjust time and sampling interval according to deployment!!!!!

%% read whole HDF file
disp('loading data...')
T = hdf5read(hdfFile, '/data/temperature');
T = T';
zVec = abs(hdf5read(hdfFile, '/dims/depth'));
    n = hdf5read(hdfFile, '/dims/time'); % raw sensor time

% by-pass, problem...
if isempty(varargin) == 1
    whichMooring = [];
    %    n = hdf5read(hdfFile, '/dims/time'); % raw sensor time
    t0_str = hdf5read(hdfFile, '/', 'StartDate');
    t0 = datenum(t0_str.Data, 'dd-mm-yyyy HH:MM:SS');
    samplingInterval = hdf5read(hdfFile, '/', 'SamplingInterval');
    timeSensor = [t0:samplingInterval/86400:t0+samplingInterval/86400*(length(n)-1)]';
    clear n 
else
    % For S1:
    whichMooring = varargin{1};
    if strcmp(whichMooring, 'S1') == 1
        disp('Working on S1, correct t0')
        t0 = datenum(2010, 2, 28, 11, 02, 0); 
        samplingInterval = 0.2;
    % For T1C    
    elseif strcmp(whichMooring, 'T1C') == 1 % For T1C:
        disp('Working on T1C, correct t0')
        t0 = datenum(2010,1,1)+62.625;
        samplingInterval = 0.2;
     % For T1B    
    elseif strcmp(whichMooring, 'T1B') == 1 % For T1B:
        disp('Working on T1B, correct t0')
        t0 = datenum(2010, 3, 7, 11, 16, 46); % <=== check this!
        t0 = datenum(2010, 3, 7, 9, 55, 45); % <=== check this! (I
                                             % fucked it up in Rostock)
        samplingInterval = 0.2;
    else
        disp('Problem! which mooring to process?')
        return
    end
    timeSensor = [t0:samplingInterval/86400:t0+samplingInterval/86400*(length(n)-1)]';
end


%% Prepare data
% make sure matrix is flip correctly
[zVec, I] = sort(zVec);
T = T(I,:);

% Extra cleaning
if strcmp(whichMooring, 'S1') == 1
    T(144:end, 1.42e6:end) = NaN;
    T(69, :) = NaN;
end

% reduce number of pts according to nskip
if nSkip ~=0
    T = T(:,1:nSkip:end);
    timeSensor = timeSensor(1:nSkip:end);
end
I = find(T<-10);
T(I) = NaN;
disp('  -> done!')



%% Nan empty values
disp('fill empty counts (NaNs from despike)...')
TFilled = hst_fillNaN(T, timeSensor);
disp('  -> done!')


%%%% Interpolate and save data %%%%
%time interpolation
disp('Time interpolation of despiked data...')
Titp = nan(size(T));
for i = 1:length(zVec)
    TVec = TFilled(i,:);    
    II = find(~isnan(TVec));
    if ~isempty(II)
        Titp(i,:) = interp1(timeSensor(II), TVec(II), timeSensor);
    end        
end
disp('  -> done!')


% vertical interpolation
disp('Vectical interpolation for missing sensors...')
Titp2 = nan(size(T));
for i = 1:length(timeSensor)
    if mod(i, 1000) == 0
        disp(sprintf('%d / %d',i,length(timeSensor)))
    end
    TVec = Titp(:,i);
    II = find(~isnan(TVec) & TVec~=0);
    Titp2(:,i) = interp1(zVec(II), TVec(II), zVec);
end
disp('  -> done!')
Titp = Titp2;

disp('Saving mat file...')
save(matFile, 'T', 'Titp', 'timeSensor', 'zVec')
disp('  -> done!')





