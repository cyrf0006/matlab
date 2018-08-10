clear all
close all

load apef_thorpe_wholeTM.mat
totalDepth = 919;


t0 = timeSensor(1);
dt = timeSensor(2)-timeSensor(1);
dt = dt*10;
rangeWindow = 1/24;
lengthWindow = round(rangeWindow./dt);
V = 5:.04:10;
tic
for i = 1:length(timeSensor)

    windowT = nan(length(zVecReg), lengthWindow);
    
    XLIM = [t0 t0+rangeWindow];
    
    I1 = find(timeSensor >= t0 & timeSensor <= t0+rangeWindow);
    
    if isempty(I1)
        windowT = nan(length(zVecReg), lengthWindow);
        windowTime = nan(1,size(windowT,2));
    else
        windowTime = timeSensor(I1);
        windowT = Titp(:,I1);
    end
    t0 = t0+dt;
    
    
    %% plot
    % *********************** Adjust_space.m ************************ %
    % Fields required by the function adjust_space.m. Please fill every
    % of the following and call "adjust_space" in the script whenever
    % you want. Do not touch four last fields
    ncol = 1; % no. subplot column
    nrow = 1; % no. subplot row
    dx = 0.03 ; % horiz. space between subplots
    dy = 0.04; % vert. space between subplots
    lefs = 0.12; % very left of figure
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

    imagesc(windowTime, totalDepth-zVecReg, windowT)
    hold on
    contour(windowTime, totalDepth-zVecReg, windowT, V, 'color', 'k')
    hold off

    %    c = colorbar('location', 'northoutside');
    %caxis([6.5 8.5])
    %ylabel(c, '\sigma_1 (kg m^{-3})', 'FontSize', 10); 
    %    ti = ylabel(c, 'T (^{\circ}C)', 'FontSize', 10, 'fontWeight','bold'); 
    ylabel('hab (m)')
    %caxis([1032 1032.1]-1000)
    datetick
    %set(gca, 'xticklabel', []);
    %set(gca, 'yticklabel', []);
    set(gca, 'ydir', 'normal');
    set(gca, 'tickdir', 'out')
    xlim([XLIM])
    ylim([0 max(totalDepth-zVecReg)])
    xlabel(sprintf(datestr(windowTime(1), 1)));
    adjust_space
    
    outfile = sprintf('window%06d.png',i);
    disp(outfile)
    print('-dpng', '-r200', outfile)
    
    close all
end

toc