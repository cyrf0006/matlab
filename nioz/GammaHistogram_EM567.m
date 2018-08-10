clear

% generated (manually) in APEF_vs_Thorpe.m, i.e.:
% save gammaInfo.mat timeSensor1 JbVec JbVec1 JbVec2 epsVec epsVec1 epsVec2 Reb
load gammaInfo2.mat

%% Choose which Gamma
GamVec = JbVec./epsVec;
%GamVec = JbVec1./epsVec1;
%GamVec = JbVec2./epsVec2;
timeVec = timeSensor1;

time2 = time2maxV('./roc12.mat',timeVec); 

% Average temperature field relative to max vel
dtAve = .2;
timeAve = [-.4:dtAve:.4];
clear GamStruct
for i = 1:length(timeAve)
    I = find(time2>=timeAve(i)-dtAve & time2<=timeAve(i)+dtAve);
    Gam = GamVec(I);
    
    %remove NaNs
    I = find(isnan(Gam));
    Gam(I) = [];
    
    varname =  sprintf('Gam%d', i);
    eval(['GamStruct.' varname '= Gam']);    
    
end
 

% Plotting section
figure(1)
clf
set(gcf,'PaperUnits','centimeters', 'paperposition', [0 0 18 14])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 5; % no. subplot column
nrow = 2; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.08; % very left of figure
rigs = 0.05; % very right of figure
tops = -0.15; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %


classEdge = [0:.05:.5];
dclass = classEdge(2)-classEdge(1);
E = struct2cell(GamStruct);
EE = [E{1} E{2} E{3} E{4} E{5}]; % whole timeserie

x = 0:pi/10:2*pi;
subplot(2,5,[1 5])
plot(x,-cos(x), 'k', 'linewidth',2)
hold on
plot(x, x*0, '--k')
set(gca, 'xticklabel', [])
set(gca, 'yticklabel', [])
ylabel('Arbitrary tidal phase')
adjust_space
pos1 = get(gca, 'pos');
adjust_space
adjust_space
adjust_space
adjust_space
pos2 = get(gca, 'pos');
pos3 = pos2;
pos3(1) = pos1(1);
pos3(3) = pos2(1)+pos2(3)-pos1(1);
pos3(4) = pos3(2)*.5;

set(gca, 'tickdir', 'out')
set(gca, 'pos', pos3);
set(gca, 'xtick', [0.2*pi:.4*pi:2*pi]);
xlim([0 2*pi])

text(pi, 0.4, 'Upslope', 'verticalAlignment', 'bottom', 'horizontalAlignment', 'center', 'fontSize', 12, 'fontWeight', 'bold')
text(pi, 0.15, '(bore/front propag.)', 'verticalAlignment', 'bottom', 'horizontalAlignment', 'center', 'fontSize', 12, 'fontWeight', 'bold')
text(2*pi-pi/100, -0.3, 'Downslope', 'verticalAlignment', 'top', 'horizontalAlignment', 'right', 'fontSize', 12, 'fontWeight', 'bold')
text(0+pi/100, -0.3, 'Downslope', 'verticalAlignment', 'top', 'horizontalAlignment', 'left', 'fontSize', 12, 'fontWeight', 'bold')

for i = 1:length(timeAve)

    subplot(2,5,i+5)
    
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
    %plot([M M], [0 .22], '--', 'color', [1 1 1]*0, 'linewidth', 1)
    %    plot([.2 .2], [0 .55], '--', 'color', [1 1 1]*0, 'linewidth', 2)
    hold off
    ylim([0 .27]);
    xlim([classEdge(1) classEdge(end)]);
    set(gca, 'xtick', classEdge(2:end-1))
    set(gca, 'ytick', [0:.1:.5])
    % text(.48, .35, sprintf('\\phi = [%1.1f,%1.1f] day', timeAve(i)-dtAve/2, timeAve(i)+dtAve/2), 'HorizontalAlignment', 'right')

    text(.48, .23, sprintf('<\\gamma> = %3.2f', M), 'HorizontalAlignment', 'right', 'fontWeight', 'bold')

    
    if i ~= 1
        set(gca, 'yticklabel', '')
    end

    xlabel('\gamma')

    set(gca, 'xtick', [0:.15:.45])
    set(gca, 'tickdir', 'out')
    %    set(gca,'xminortick','on');
    
    if i == 1
        ylabel('Relative Freq.')
    end
    
    adjust_space
end

% $$$ subplot(3,2,6)
% $$$ set(gca, 'visible', 'off')
% $$$ text(.9, .9, 'Whole timeseries stats', 'HorizontalAlignment', 'right')
% $$$ text(.9, .8, '(8-16 Oct. 2012):', 'HorizontalAlignment', 'right')
% $$$ 
% $$$ text(.9, .65, sprintf('<\\gamma> = %3.2f', nanmean(GamVec)), 'HorizontalAlignment', 'right')
% $$$ %text(.9, .525, sprintf('\\nu_{1/2} = %3.2f', median(EE)),'HorizontalAlignment', 'right')
% $$$ text(.9, .5, sprintf('\\mu_{1/2} = %3.2f', median(GamVec)),'HorizontalAlignment', 'right')
% $$$ text(.9, .3, sprintf('n = %2.0f', length(GamVec)),'HorizontalAlignment', 'right')


print('-dpng', '-r300', 'GammaHistogram_EM.png')

