function hst_TRrel(ctdList, timeFile, TLIM)
    
% function hst_TSrel(ctdList, timeFile)
%
% Made specificly to read and process CTD casts from ROC12
% cruise. Shell preccessing were done to clean file names and to
% build the timeFile
%    
% usage ex:
% hst_TRrel('ctdFiles.list','./CTD360/fred_processed/timeFile.txt', [6 9])
% hst_TRrel('ctdFiles24h.list','./CTD360/fred_processed/timeFile24h.txt', [6 9])
% hst_TRrel('ctdFilesMooring.list','./CTD360/fred_processed/timeFileMooring.txt', [6.5 9]) <-----GOOD ONE




% load CTD files
fid = fopen(ctdList);
C = textscan(fid, '%s', 'delimiter', '\n');
ctd = char(C{1});
noFiles = size(ctd,1);

% get time vector
fid = fopen(timeFile);
C = textscan(fid, '%s', 'delimiter', '\n');
ctdTime = char(C{1});

TVec = [];
ZVec = [];
SVec = [];
RVec = [];
timeVec = [];

for i = 1:noFiles

    fname = ctd(i, :);
    I = find(fname==' ');   
    fname(I) = [];
    
    data = load(fname);
   
% $$$     %% old version
% $$$ % $$$     Z = data(:, 14);
% $$$ % $$$     T = data(:, 15);
% $$$ % $$$     S = data(:, 16);'./CTD360/fred_processed/
% $$$     Z = data(:, 1);
% $$$     T = data(:, 2);
% $$$     S = data(:, 15);
% $$$     %    S = gsw_SP_from_C(C,T,Z);
% $$$ 
% $$$      % Isolate downcast
% $$$     dz = gradient(Z);
% $$$     I = find(dz>0.5*mean(abs(dz)));
% $$$     
% $$$     %T,S -> CT,SA & rho
% $$$     SA = gsw_SA_from_SP(S(I), Z(I), -15.1, 55.5);
% $$$     CT = gsw_CT_from_t(SA,T(I),Z(I));
% $$$     
% $$$     rho = gsw_rho(SA,CT,Z(I));     
% $$$ 
% $$$     ZVec = [ZVec; Z(I)];
% $$$     TVec = [TVec; CT];
% $$$     SVec = [SVec; SA];
% $$$     RVec = [RVec;rho];        
% $$$     timeVec = [timeVec; Z(I)*0+ctdTime(i)];    
    
    
    %% new version
    Z = data(:, 1);
    %    T = data(:, 2);
    T = data(:, 15); %(potential)
    rho = data(:, 11);    

    % Isolate downcast
    dz = gradient(Z);
    I = find(dz>0.5*mean(abs(dz)));
    
    ZVec = [ZVec; Z(I)];
    TVec = [TVec; T(I)];
    SVec = [SVec; rho(I)];
    RVec = [RVec;rho(I)];        
    timeVec = [timeVec; Z(I)*0+ctdTime(i)];
    
end


% Clean T,S
I = find(SVec<1e-6);
SVec(I) = [];
TVec(I) = [];
RVec(I) = [];
ZVec(I) = [];
I = find(TVec<1e-6);
SVec(I) = [];
TVec(I) = [];
RVec(I) = [];
ZVec(I) = [];
I = find(RVec<1e-6);
SVec(I) = [];
TVec(I) = [];
RVec(I) = [];
ZVec(I) = [];


figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 10])
I = find(TVec>=TLIM(1) & TVec <= TLIM(2));
TVec = TVec(I);
RVec = RVec(I);

[Y, I] = sort(TVec);
YY  = RVec(I);
p = polyfit(TVec(I), RVec(I), 1);
Tfit = min(TVec):.1:max(TVec);
Rfit = Tfit*p(1) + p(2);

plot(TVec, RVec, '.k')
hold on
plot(Tfit, Rfit, 'color', [.5 .5 .5]) 

Rfit = Y*p(1) + p(2);
R = corrcoef(Rfit, RVec(I));
title(sprintf('[%2.1f-%2.1f]degC / R=%3.2f', TLIM(1), TLIM(2), R(2)))
xlabel('\theta_0(^{\circ}C)')
ylabel('\sigma_{1} (kg m^{-3})')
print('-depsc2', 'temp_rho.eps')

save temp_rho_relFile.mat p
disp('temp_rho_relFile.mat saved!')
p(1)



%% bootstrap for error
disp('bootstrap...')
nboot = 1000;

N = length(TVec);
p_boot_b0 = nan(1,nboot);
p_boot_b1 = nan(1,nboot);
p_boot_b2 = nan(1,nboot);
p_boot_b3 = nan(1,nboot);

for i = 1:nboot
   
    r = rand(N,1);
    r = ceil(r*N/1);
    
    myVecT = TVec(r);
    myVecR = RVec(r);
    
    [Y, I] = sort(myVecT);
    p = polyfit(myVecT(I), myVecR(I), 3);
    p_boot_b3(i) = p(1);
    p_boot_b2(i) = p(2);
    p_boot_b1(i) = p(3);
    p_boot_b0(i) = p(4);
end

% $$$ p_boot_dot = nanmean(p_boot_b);
% $$$ perror = sqrt(sum([p_boot_b-p_boot_dot].^2)./(nboot-1));
% $$$ 
% $$$ [p_boot_dot  perror]

CI_2p5 = round(2.5/100*nboot);
CI_97p5 = round(97.5/100*nboot);

disp('    ave      2.5%      97.5%')
psort = sort(p_boot_b3);
[nanmean(psort)  psort(CI_2p5) psort(CI_97p5)]
psort = sort(p_boot_b2);
[nanmean(psort)  psort(CI_2p5) psort(CI_97p5)]
psort = sort(p_boot_b1);
[nanmean(psort)  psort(CI_2p5) psort(CI_97p5)]
psort = sort(p_boot_b0);
[nanmean(psort)  psort(CI_2p5) psort(CI_97p5)]

% standard deviation of error
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 2; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.08 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.1; % very right of figure
tops = 0.05; % top of figure
bots = 0.18; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

figure(2)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 16])
subplot(121)
I = find(TVec>=TLIM(1) & TVec <= TLIM(2));
TVec = TVec(I);
RVec = RVec(I);

[Y, I] = sort(TVec);
YY  = RVec(I);
p = polyfit(TVec(I), RVec(I), 3);
Tfit = min(TVec):.1:max(TVec);
%Rfit = Tfit*p(1) + p(2);
Rfit = polyval(p, Y);

plot(RVec, TVec, '.k')
hold on
%plot(Rfit, Tfit, 'color', [.5 .5 .5]) 
plot(Rfit, Y, 'color', [.5 .5 .5], 'lineWidth', 2) 

%Rfit = Y*p(1) + p(2);
R = corrcoef(Rfit, RVec(I));
%title(sprintf('[%2.1f-%2.1f]degC / R=%3.2f', TLIM(1), TLIM(2), R(2)))
set(gca, 'tickdir', 'out')
ylabel('\theta_0(^{\circ}C)', 'fontSize', 14, 'fontWeight', 'bold')
XLAB = xlabel('\sigma_{1} (kg m^{-3})', 'fontSize', 14, 'fontWeight', 'bold');
ylim(TLIM)
xlim([31.8 32.2])
text(31.85, 6.6, 'a', 'fontSize', 16, 'fontWeight', 'bold')
%text(31.95, 8.75, '\sigma_{fit} = a_3T^3 + a_2T^2 + a_1T^1 + a_0', 'fontSize', 9)
adjust_space
XLAB_POS = get(XLAB, 'pos');
XLAB_POS(2) = 6.2;
set(XLAB, 'pos', XLAB_POS);


subplot(122)
plot(YY-Rfit, Y, '.k')
hold on
plot([0 0], [min(Y) max(Y)], 'color', [.5 .5 .5], 'lineWidth', 2)
%xlabel('$\rm \sigma_{1000}-\,~\tilde{\sigma}_{1000}~(kg~m^{-3})$', 'interpreter', 'latex')
XLAB = xlabel('\sigma_{1_{obs}}-{\sigma}_{1_{fit}} (kg m^{-3})', 'fontSize', 14, 'fontWeight', 'bold');
set(gca, 'tickdir', 'out')
set(gca, 'yticklabel', [])
ylim(TLIM)
xlim([-.02 .025])
text(-.015, 6.6, 'b', 'fontSize', 16, 'fontWeight', 'bold')
adjust_space
XLAB_POS = get(XLAB, 'pos');
XLAB_POS(2) = 6.2;
set(XLAB, 'pos', XLAB_POS);

set(gcf, 'renderer', 'painters')
print('-depsc2', 'temp_rho_3rd.eps')
save temp_rho_relFile_3rd.mat p
disp('temp_rho_relFile_3rd.mat saved!')
