function ctd_timeseries(ctdList, timeFile)
    
% function hst_TSrel(ctdList, timeFile)
%
% Made specificly to read and process CTD casts from ROC12
% cruise. Shell preccessing were done to clean file names and to
% build the timeFile. I also removed few lines in one CTD file
% (64PE360_S101C01fil_ctm_drv_bav_buo.asc)
    
%    
% usage ex: ctd_timeseries('ctdFilesYoYo.list','./CTD360/fred_processed/timeFileYoYo.txt')
%  if run in /home/cyrf0006/research/NIOZ/RockallBank
g = 9.81;
zMax = 796.7500;
zMin = 915.7500;

% load CTD files
fid = fopen(ctdList);
C = textscan(fid, '%s', 'delimiter', '\n');
ctd = char(C{1});
noFiles = size(ctd,1);

% get time vector
fid = fopen(timeFile);
C = textscan(fid, '%s', 'delimiter', '\n');
ctdTime = char(C{1});
ctdTime = datenum(ctdTime);

TVec = [];
ZVec = [];
SVec = [];
N2Vec = [];
O2Vec = [];
RVec = [];
timeVec = [];
for i = 1:noFiles

    fname = ctd(i, :);
    I = find(fname==' ');   
    fname(I) = [];
    
    data = load(fname);
   
    Z = data(:, 1);
    T = data(:, 2);
    C = data(:, 3);
    O2 = data(:, 6);

    SP = gsw_SP_from_C(C,T,Z);
    SA = gsw_SA_from_SP(SP, Z,-13.8, 57.5);
    CT = gsw_CT_from_t(SA, T, Z);
    sig1 = gsw_sigma1(SA,CT);
    rho = gsw_rho(SA,CT,Z);

    % Isolate downcast
    dz = gradient(Z);
    I = find(dz>0.5*mean(abs(dz)));    
% $$$     [zm, I] = max(Z);
% $$$     I = 1:I;
    
    rhoSort = sort(rho(I));
    N2 = g./nanmean(rho(I)).*gradient(rhoSort)./gradient(Z(I));
    
    ZVec = [ZVec; NaN; Z(I)];
    TVec = [TVec; NaN; CT(I)];
    RVec = [RVec; NaN; sig1(I)];        
    N2Vec = [N2Vec; NaN; N2]; 
    O2Vec = [O2Vec; NaN; O2];
    timeVec = [timeVec; NaN; Z(I)*0+ctdTime(i)];        
end

% grid the data
myZ = 1:10:round(max(ZVec));
myTime = min(timeVec):.02:max(timeVec);
[X, Y] = meshgrid(myTime, myZ);
I = find(~isnan(TVec));
[XI,YI,Tgrid] = griddata(timeVec(I), ZVec(I), TVec(I),X,Y);
[XI,YI,Rgrid] = griddata(timeVec(I), ZVec(I), RVec(I),X,Y);
[XI,YI,O2grid] = griddata(timeVec(I), ZVec(I), O2Vec(I),X,Y);

I = find(~isnan(N2Vec));
[XI,YI,N2grid] = griddata(timeVec(I), ZVec(I), N2Vec(I),X,Y);


figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 15 18])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 2; % no. subplot row
dx = 0.0; % horiz. space between subplots
dy = 0.02; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.14; % very right of figure
tops = 0.03; % top of figure
bots = 0.07; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %


subplot(211)
contourf(XI, YI, Tgrid, 100, 'linestyle', 'none')
set(gca, 'ydir', 'reverse')
datetick
xlim([min(myTime) max(myTime)])
ylabel('Depth (m)')
set(gca, 'xticklabel', [])
set(gca, 'tickdir','out')
hold on
plot(timeVec, ZVec, '--k')
plot([min(timeVec) max(timeVec) max(timeVec) min(timeVec) min(timeVec)], ...
     [zMax zMax zMin zMin zMax], '--k', 'linewidth', 2)
hold off
cb = colorbar;
adjust_space
drawnow 

cbPos = get(cb, 'pos');
cbPos(1) = cbPos(1)-dx;
cbPos(2) = cbPos(2)+dy;
cbPos(4) = cbPos(4)-1*dy;
cbPos(3) = cbPos(3)*.6;
set(cb, 'pos', cbPos)
ti = ylabel(cb,'T (^{\circ}C)', 'FontSize', 10);
drawnow


subplot(212)
contourf(XI, YI, log10(N2grid), 20, 'linestyle', 'none')
%imagesc(XI(1,:), YI(:,1), log10(N2grid)); %shading interp;
set(gca, 'ydir', 'reverse')
set(gca, 'tickdir','out')
datetick
xlim([min(myTime) max(myTime)])
xlabel(sprintf('%s,%s', datestr(myTime(1), 7), datestr(myTime(end), 1)))
hold on
contour(XI, YI, log10(Rgrid), 20, 'k')
plot(timeVec, ZVec, '--k')
plot([min(timeVec) max(timeVec) max(timeVec) min(timeVec) min(timeVec)], ...
     [zMax zMax zMin zMin zMax], '--k', 'linewidth', 2)
hold off
cb = colorbar;
caxis([-5.5 -4]) 
adjust_space
drawnow 


cbPos = get(cb, 'pos');
cbPos(1) = cbPos(1)-dx;
cbPos(2) = cbPos(2)+dy;
cbPos(4) = cbPos(4)-1*dy;
cbPos(3) = cbPos(3)*.6;
set(cb, 'pos', cbPos)
ti = ylabel(cb,'log(N^2 / (s^{1}))', 'FontSize', 10);
drawnow

print('-dpng', '-r300', './figures/ctd_timeserie.png')
set(gcf, 'renderer', 'painters')
print('-depsc2', './ctd_timeserie.eps')
keyboard