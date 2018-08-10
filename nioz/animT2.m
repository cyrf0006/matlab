clear all
close all

load apef_thorpe_wholeTM.mat
totalDepth = 919;


%% play with ADCP vel.
load('./roc12.mat');

N = [SerNmmpersec/1000]'; % in m/s now
E = [SerEmmpersec/1000]';
zAdcp = AnDepthmm/1000; % instrument depth
timeAdcp = [datenum(SerYear+2000, SerMon, SerDay, SerHour, SerMin, SerSec)]';

% isolate good timeserie (remove if depth <10)
I = abs(zAdcp - nanmean(zAdcp))>10;
zAdcp(I) = [];
timeAdcp(I) = [];
N(:,I) = [];
E(:,I) = [];
clear Ser* An*

% eliminate data below bottom
I = find(abs(E)>10);
E(I) = NaN;
I = sum(isnan(E),2);
II = find(I>mean(I));
Iremove = II(1):length(I);
E(Iremove,:) = [];
N(Iremove,:) = [];
zAdcp(Iremove) = [];

[U,V] = rotate_vecd(E,N,-30);
uVec = nanmean(U,1);
vVec = nanmean(V,1);

% Filter timeserie
dt = diff(abs(timeAdcp(1:2))); %day
fs = 1/dt; % cpd
freq_low = 1; %cpd
Wn_low = freq_low/(fs/2);
[b,a] = butter(4, Wn_low);
Vfilt = filtfilt(b,a,vVec);

II = find(timeAdcp>=datenum(2012,10,10,0,0,0) & timeAdcp<=datenum(2012,10,11,23,0,0));
Vfilt = Vfilt(II);
timeAdcp = timeAdcp(II);

%% Now with T data
I = find(timeSensor>=datenum(2012,10,10,20,0,0) & timeSensor<=datenum(2012,10,11,5,0,0));
Titp = Titp(:,I);
timeSensor = timeSensor(I);

t0 = timeSensor(1);
dt = timeSensor(2)-timeSensor(1);
dt = dt*45; % time step between frames
rangeWindow = 1/24; % length of the window
lengthWindow = round(rangeWindow./dt);
V = 5:.02:10;
tic
for i = 1:length(timeSensor)/45
     %% plot
    % *********************** Adjust_space.m ************************ %
    % Fields required by the function adjust_space.m. Please fill every
    % of the following and call "adjust_space" in the script whenever
    % you want. Do not touch four last fields
    ncol = 1; % no. subplot column
    nrow = 4; % no. subplot row
    dx = 0.03 ; % horiz. space between subplots
    dy = 0.04; % vert. space between subplots
    lefs = 0.14; % very left of figure
    rigs = 0.05; % very right of figure
    tops = 0.03; % top of figure
    bots = 0.1; % bottom of figure
    figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
    figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
    count_col = 1;
    count_row = 1;
    % *************************************************************** %
    figure('visible', 'off')
    clf
    set(gcf, 'PaperUnits', 'Centimeters', 'PaperPosition', [6 1 12 12])
    
    I1 = find(timeSensor >= t0 & timeSensor <= t0+rangeWindow);
    
    windowT = nan(length(zVecReg), lengthWindow);
    
    XLIM = [t0 t0+rangeWindow];
    
    if isempty(I1)
        windowT = nan(length(zVecReg), lengthWindow);
        windowTime = nan(1,size(windowT,2));
    else
        windowTime = timeSensor(I1);
        windowT = Titp(:,I1);
    end
    t0 = t0+dt;
    
    %% S1
    subplot(411)
    Vd = detrend(Vfilt);
    plot(timeAdcp, Vd, 'k', 'linewidth', 2)
    hold on
    plot(timeAdcp, Vd*0, '--k', 'linewidth', 1)
    xlim([datenum(2012,10,10,17,0,0) datenum(2012,10,11,9,0,0)])
    datetick('x', 15, 'keepticks')
    %    xlim([datenum(2012,10,10,6,0,0) datenum(2012,10,11,6,0,0)])
    hold on
    [Y, I] = min(abs(mean(timeSensor(I1))-timeAdcp));
    plot(timeAdcp(I), Vd(I), '.r', 'markersize', 25)
    hold off
    ylabel('V_{filt} (m/s)')
    adjust_space
    
    %% S2
    subplot(4,1,[2 4])
    imagesc(windowTime, totalDepth-zVecReg, windowT)
    hold on
    contour(windowTime, totalDepth-zVecReg, windowT, V, 'color', 'k')
    hold off

    %    c = colorbar;
    caxis([6.8 8.5])
    %ylabel(c, '\sigma_1 (kg m^{-3})', 'FontSize', 10); 
    %    ti = ylabel(c, 'T (^{\circ}C)', 'FontSize', 10, 'fontWeight','bold'); 
    ylabel('hab (m)')
    %caxis([1032 1032.1]-1000)    
    set(gca, 'xticklabel', []);
    %set(gca, 'yticklabel', []);
    set(gca, 'ydir', 'normal');
    set(gca, 'tickdir', 'out')
    xlim([XLIM])
    ylim([0 max(totalDepth-zVecReg)])
    %xlabel(sprintf(datestr(windowTime(1), 1)));
    xlabel('10-11 Oct. 2012')
    adjust_space
    
    pos1 = get(gca, 'pos');
    adjust_space
    adjust_space
    pos2 = get(gca, 'pos');
    pos2(4) = 3*pos2(4)+dy;
    set(gca, 'pos', pos2)
    
    
    outfile = sprintf('window%06d.png',i);
    disp(outfile)
    print('-dpng', '-r300', outfile)
    close all
end

toc