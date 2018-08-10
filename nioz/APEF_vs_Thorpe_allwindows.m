clear all
close all

load apef_thorpe_wholeTM.mat
totalDepth = 919;
   
Ri = hst_Ri('./roc12.mat',N2bg,zVecReg, timeSensor);

for ii = 17:227
%for ii = 164

    I1 = find(timeSensor >= datenum(2012,10,7)+ii/(24) & timeSensor <= datenum(2012,10,7,1,0,0)+ii/(24));  
    I2 = find(Ri.timeVec >= datenum(2012,10,7)+ii/(24) & Ri.timeVec <= datenum(2012,10,7,1,0,0)+ii/(24));  

    if isempty(I1)
        continue
    end
    
    windowRho = Ritp1(:,I1);
    windowTime = timeSensor(I1);
    windowT = Titp(:,I1);

    H = size(windowRho,1);

    Re = hst_Re('./roc12.mat',zVecReg, windowTime);
    
    JbVec = []; 
    epsVec = [];
    for j = 1:length(windowTime)
        rhoVec_raw =  windowRho(:,j);
        [rhoVec_sort, I] = sort(rhoVec_raw);
        d = zVecReg - zVecReg(I);
        EP1 = sum(rhoVec_raw.*(max(zVecReg)-zVecReg')); 
        EP2 = sum(rhoVec_sort.*(max(zVecReg)-zVecReg'));
        rho0 = nanmean(windowRho(:,j));

        % local stratification
        N2local = g./rho0.*gradient(rhoVec_sort)./gradient(zVecReg');
        JbVec = [JbVec g./H/rho0.*(EP1-EP2).*nanmean(sqrt(N2local))];
        epsVec = [epsVec .64*rms(d).^2.*nanmean(sqrt(N2local)).^3];
    end

    
    %% Few plots.
    % *********************** Adjust_space.m ************************ %
    % Fields required by the function adjust_space.m. Please fill every
    % of the following and call "adjust_space" in the script whenever
    % you want. Do not touch four last fields
    ncol = 1; % no. subplot column
    nrow = 4; % no. subplot row
    dx = 0.03 ; % horiz. space between subplots
    dy = 0.04; % vert. space between subplots
    lefs = 0.23; % very left of figure
    rigs = 0.14; % very right of figure
    tops = 0.06; % top of figure
    bots = 0.07; % bottom of figure
    figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
    figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
    count_col = 1;
    count_row = 1;
    % *************************************************************** %
    figure(1)
    clf
    set(gcf, 'PaperUnits', 'Centimeters', 'PaperPosition', [6 1 16 ...
                        18], 'paperOrientation', 'landscape')
    s1 = subplot(4,1,[1 2]);
    %imagesc(windowTime, zVecReg, Ritp1-1000)
    imagesc(windowTime, totalDepth-zVecReg, windowT)
    hold on
    contour(windowTime, totalDepth-zVecReg, windowT, 20, 'color', 'k')
    hold off
    c = colorbar('location', 'northoutside');
    %ylabel(c, '\sigma_1 (kg m^{-3})', 'FontSize', 10); 
    ti = ylabel(c, 'T (^{\circ}C)', 'FontSize', 10, 'fontWeight','bold'); 
    %ylabel('hab (m)')
    %caxis([1032 1032.1]-1000)
    datetick
    set(gca, 'xticklabel', []);
    set(gca, 'yticklabel', []);
    set(gca, 'ydir', 'normal');
    set(gca, 'tickdir', 'out')
    xlim([windowTime(1) windowTime(end)])
    ylim([0 max(totalDepth-zVecReg)])

    text(windowTime(1)+.5/1440, 10, sprintf('Re = %2.1e', nanmean(Re)), 'horizontalAlignment', 'left','verticalAlignment','bottom','color',[1 1 1]*1,'fontweight','bold')

    adjust_space
    pos1 = get(gca, 'pos');
    adjust_space
    pos2 = get(gca, 'pos');
    pos2(4) = 2*pos2(4)+dy;
    set(gca, 'pos', pos2)
    drawnow

    cbPos = get(c, 'pos');
    cbPos(1) = cbPos(1)+.08;
    cbPos(3) = cbPos(3)-.1;
    cbPos(2) = cbPos(2)-.04;
    cbPos(4) = cbPos(4)*.4;
    set(c, 'pos', cbPos)
    set(ti, 'rotation', 0)
    CLIM = get(gca, 'clim');
    ti_pos = [CLIM(1)-.02, 1 1];
    set(ti, 'pos', ti_pos)
    drawnow

    s2 = subplot(4,1,3);
    semilogy(windowTime, JbVec, 'k') 
    hold on
    plot([timeSensor(1) timeSensor(end)], [1 1]*nanmean(JbVec), '--k')
    semilogy(windowTime, epsVec, 'color', [1 1 1]*.5)
    plot([timeSensor(1) timeSensor(end)], [1 1]*nanmean(epsVec), 'color',  [1 1 1]*.5, 'linestyle', '--')
    datetick
    set(gca, 'xticklabel', []);
    set(gca, 'tickdir', 'out')
    xlim([windowTime(1) windowTime(end)])
    %xlabel(sprintf(datestr(timeSensor(I(1)), 1)));
    ylabel('J_b,\epsilon (m^2 s^{-3})')
    set(gca, 'ytick', [1e-9 1e-8 1e-7 1e-6])
    ylim([1e-9 1e-5])
    adjust_space
    drawnow

    s3 = subplot(4,1,4);
    plot(windowTime,JbVec./epsVec, 'k')
    %datetick('x', 7)
    hold on
    plot(windowTime,JbVec*0+.2, 'color', [1 1 1]*.5,'linestyle', '--')
    ylim([0 1])
    datetick('x', 15)
    xlim([windowTime(1) windowTime(end)])
    xlabel(sprintf(datestr(windowTime(1), 1)));
    %xlabel('2012-10-11')
    %xlabel('October 2012')
    ylabel('\Gamma');
    text(windowTime(end)-12/1440, .75, sprintf('<\\Gamma> = %1.2f',nanmean(JbVec./epsVec)))
    set(gca, 'tickdir', 'out')
    adjust_space


    %% add Richardson number
    POS = get(s1, 'pos');
    ax = axes('pos', [0.08 POS(2) 0.12 POS(4)], 'yticklabel', [], 'xticklabel',[]);
    plot(nanmean(Ri.Ri(:,I2),2), totalDepth-Ri.zVec, 'k', 'linewidth', 2)
    hold on
    plot(Ri.zVec*0+.25, totalDepth-Ri.zVec, '--k')
    set(ax, 'ydir', 'normal');
    set(ax, 'xscale', 'log');
    set(ax, 'tickdir', 'out')
    ylim([0 max(totalDepth-zVecReg)])
    ylabel('hab (m)')
    xlabel('Ri')

    set(gcf, 'renderer', 'painters')
    outfile = sprintf('window%04d.pdf',ii);
    disp(outfile)
    print('-dpdf', outfile)

end