function time2maxV = extract_tide(adcpFile, timeVec)

% function L = extract_tide(adcpFile, zVec, timeVec)
%
% usage ex: 
%  L = extract_tide('/media/Seagate1TB/NIOZ/thermistordata/ROC12/roc12.mat')   
    
%% Load and get the relevant variables
load(adcpFile);

N = [SerNmmpersec/1000]'; % in m/s now
E = [SerEmmpersec/1000]';
W = [SerVmmpersec/1000]';

zAdcp = AnDepthmm/1000; % instrument depth
timeAdcp = [datenum(SerYear+2000, SerMon, SerDay, SerHour, SerMin, SerSec)]';

% isolate good timeserie (remove if depth <10)
I = abs(zAdcp - nanmean(zAdcp))>10;
zAdcp(I) = [];
timeAdcp(I) = [];
N(:,I) = [];
E(:,I) = [];
W(:,I) = [];

clear Ser* An*

% eliminate data below bottom
I = find(abs(E)>10);
E(I) = NaN;
I = sum(isnan(E),2);
II = find(I>mean(I));
Iremove = II(1):length(I);
E(Iremove,:) = [];
N(Iremove,:) = [];
W(Iremove,:) = [];
zAdcp(Iremove) = [];


% roation may be not necessary
[U,V] = rotate_vecd(E,N,30);
uVec = nanmean(U,1);
vVec = nanmean(V,1);

% Filter timeserie
dt = diff(abs(timeAdcp(1:2))); %day
fs = 1/dt; % cpd
freq_low = 1; %cpd
Wn_low = freq_low/(fs/2);
[b,a] = butter(4, Wn_low);
Ufilt = filtfilt(b, a, uVec);
Vfilt = filtfilt(b,a,vVec);

Vitp = interp1(timeAdcp, Vfilt, timeVec);

% find max velocities
keyboard

timeMax = [];
timeMin = [];
for i = 2:length(timeVec)
    if Vitp(i-1)<Vitp(i) & Vitp(i+1)<Vitp(i)
        timeMax = [timeMax timeVec(i)];        
    end
    
    if Vitp(i-1)>Vitp(i) & Vitp(i+1)>Vitp(i)
        timeMin = [timeMin timeVec(i)];        
    end
end

time2maxV = nan(size(timeVec));
for i = 1:length(timeVec)
    [Y,I] = min(abs(timeVec(i)-timeMax));
    time2maxV(i) = timeVec(i)-timeMax(I);
end
