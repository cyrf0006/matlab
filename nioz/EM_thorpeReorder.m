clear all
close all

load apef_thorpe_wholeTM.mat
totalDepth = 919;


% subsets:
%I1 = find(timeSensor >= datenum(2012,10,8,12,25,0)-.5 & timeSensor<= datenum(2012,10,8,13,25,0)-.5);%ovt
I1 = find(timeSensor >= datenum(2012,10,10,23,0,0) & timeSensor <= datenum(2012,10,11,0,0,0));%cvt

T1 = T(:,I1);
timeSensor = timeSensor(I1);
% Vertical interpolation for missing sensors and regular grid
Titp = [];
for i = 1:length(timeSensor)
    TVec = T1(:,i);
    II = find(~isnan(TVec) & TVec~=0);
    Titp = [Titp; [interp1(zVec(II), TVec(II), zVec)]'];
end
Titp = Titp';

%[Y, I] = min(abs(timeSensor - datenum(2012,10,8,0,59,0)));
[Y, I] = min(abs(timeSensor - datenum(2012,10,10,23,35,0)));

figure(1)
clf
set(gcf, 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 15 10])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 3; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.06 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.06; % very left of figure
rigs = 0.04; % very right of figure
tops = 0.03; % top of figure
bots = 0.11; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %
V1 = 5:.01:10;
V = 5:.05:10;
s1 = subplot(1,3,[1 2]);
contourf(timeSensor, totalDepth-zVec, Titp, V1, 'linestyle', 'none')
hold on
contour(timeSensor, totalDepth-zVec, Titp, V, 'linecolor', 'k')
plot([timeSensor(I) timeSensor(I)], totalDepth-[zVec(1) zVec(end)], '--m', 'linewidth', 2)
datetick
xlim([timeSensor(1), timeSensor(end)])
%ylim([min(zVec) max(zVec)])
ylim([0 max(totalDepth-zVecReg)])
xlabel(sprintf(datestr(timeSensor(I), 1)));
set(gca, 'ydir', 'normal')
ylabel('hab (m)')

adjust_space
POS1 = get(gca, 'pos');
adjust_space
POS2 = get(gca, 'pos');
POS3 = POS2;
POS3(1) = POS1(2);
POS3(3) = POS2(1)+POS2(3)-POS1(1);
set(gca, 'pos', POS3)


s2 = subplot(133);
plot(Titp(:,I), totalDepth-zVec, 'k')
hold on
plot(flipud(sort(Titp(:,I))), totalDepth-zVec, 'r')
set(gca, 'ydir', 'normal')
set(gca, 'yticklabel', [])
set(gca, 'fontSize', 10)
ylim([0 max(totalDepth-zVecReg)])
%ylim([min(zVec) max(zVec)])
xlabel('\theta_1(^{\circ}C)')
xlim([7 7.8])
adjust_space
print('-dpng', '-r300', 'EM_Thorpe_toberename.png')


figure(2)
clf
set(gcf, 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 15 10])
contourf(timeSensor, totalDepth-zVec, Titp, V1, 'linestyle', 'none')
hold on
contour(timeSensor, totalDepth-zVec, Titp, V, 'linecolor', 'k')
hold off
set(gca, 'ydir', 'normal')
datetick
xlim([timeSensor(1), timeSensor(end)])
ylim([0 max(totalDepth-zVecReg)])
%ylim([min(zVec) max(zVec)])
xlabel(sprintf(datestr(timeSensor(I), 1)));
ylabel('hab (m)')
xlabel(sprintf(datestr(timeSensor(I), 1)));
print('-dpng', '-r300', 'EM_temp_toberename.png')
