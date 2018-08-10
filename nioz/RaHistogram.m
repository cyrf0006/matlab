clear

% generated (manually) in APEF_vs_Thorpe.m, i.e.:
% save gammaInfo.mat timeSensor1 JbVec JbVec1 JbVec2 epsVec epsVec1 epsVec2 Reb
load gammaInfo2.mat

%% Choose which Gamma
%GamVec = JbVec./epsVec;
%GamVec = JbVec1./epsVec1;
%GamVec = JbVec2./epsVec2;
timeVec = timeSensor1;
time2 = time2maxV('./roc12.mat',timeVec); 

allAve = nanmean(log10(RaVec2));

% Average temperature field relative to max vel
dtAve = .2;
timeAve = [-.4:dtAve:.4];
clear RaStruct
for i = 1:length(timeAve)
    I = find(time2>=timeAve(i)-dtAve & time2<=timeAve(i)+dtAve);
    Ra = log10(RaVec2(I));
    
    %remove NaNs
    I = find(isnan(Ra));
    Ra(I) = [];
    
    varname =  sprintf('Ra%d', i);
    eval(['RaStruct.' varname '= Ra']);    
    
end
 

% Plotting section
figure(1)
clf
set(gcf,'PaperUnits','centimeters', 'paperposition', [0 0 14 18])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 2; % no. subplot column
nrow = 3; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.03; % vert. space between subplots
lefs = 0.12; % very left of figure
rigs = 0.05; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %


classEdge = [4:.5:16];
dclass = classEdge(2)-classEdge(1);
E = struct2cell(RaStruct);
EE = [E{1} E{2} E{3} E{4} E{5}]; % whole timeserie

for i = 1:length(timeAve)

    subplot(3,2,i)
    
    EVec = E{i};
    I = find(EVec>max(classEdge));
    EVec(I) = max(classEdge);
    
    clear relativeNo;
    for j = 1:length(classEdge)-1
        I = find(EVec>=classEdge(j) & EVec<=classEdge(j+1));
        relativeNo(j) = length(I)/length(EVec);
    end
  
    bar(classEdge(1:end-1)+dclass/2, relativeNo, 'facecolor', [.2 .2 .2]);
    hold on
    M = nanmean(E{i});
    plot([M M], [0 .22], '--', 'color', [1 1 1]*0, 'linewidth', 1)
    plot([allAve allAve], [0 .22], '--', 'color', [1 1 1]*0, 'linewidth', 2)
    hold off
    ylim([0 .2]);
    xlim([classEdge(1) classEdge(end)]);
    set(gca, 'xtick', classEdge(2:end-1))
    set(gca, 'ytick', [0:.1:.5])
    text(15.9, .16, sprintf('\\phi = [%1.1f,%1.1f] day', timeAve(i)-dtAve/2, timeAve(i)+dtAve/2), 'HorizontalAlignment', 'right')

    text(15.9, .13, sprintf('<log_{10}(Ra)> = %3.2f', M), 'HorizontalAlignment', 'right')
% $$$     str = sprintf('$$\\tilde{\\Gamma} = %3.2f$$', median(E{i}));
% $$$     text('Interpreter','latex','Position',[.35 .25],'String',str, 'HorizontalAlignment', 'right')
    text(15.9, .10, sprintf('\\mu_{1/2} = %3.2f', median(E{i})),'HorizontalAlignment', 'right')

    
    if i == 2 | i == 4 | i == 6
        set(gca, 'yticklabel', '')
    end
    if i == 5 
        xlabel('log_{10}(Ra)')
    end
    set(gca, 'xtick', [4:2:16])
    set(gca, 'tickdir', 'out')
    %    set(gca,'xminortick','on');
    
    if i == 1 | i == 2 | i == 3 %| i == 4
        set(gca, 'xticklabel', '')
    else
        %        set(gca, 'xticklabel', ['0', '', '0.1', '', '0.2', '', '0.3', ...
        %          '', '0.4', '', '0.5'])
        %        minorTicks = [0:.05:.5];
        %        minorTicks(2:2:end) = sparse;
        
% $$$         set(gca, 'xticklabel',minorTicks)
    end
    if i == 3
        ylabel('Relative Freq.')
    end
    
    adjust_space
end

subplot(3,2,6)
set(gca, 'visible', 'off')
text(.9, .9, 'Whole timeseries stats', 'HorizontalAlignment', 'right')
text(.9, .8, '(8-16 Oct. 2012):', 'HorizontalAlignment', 'right')

text(.9, .65, sprintf('<log_{10}(Ra)> = %3.2f', nanmean(log10(RaVec2))), 'HorizontalAlignment', 'right')
%text(.9, .525, sprintf('\\nu_{1/2} = %3.2f', median(EE)),'HorizontalAlignment', 'right')
text(.9, .5, sprintf('\\mu_{1/2} = %3.2f', median(log10(RaVec2))),'HorizontalAlignment', 'right')
text(.9, .3, sprintf('n = %2.0f', length(RaVec2)),'HorizontalAlignment', 'right')


print('-depsc2', 'RaHistogram_tidal.eps')

