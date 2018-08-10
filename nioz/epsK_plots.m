clear all
close all
FS = 12;
FS2 = 14;
nu_cst = 1e-6; %m^2/s
kappa_T = 1.4e-7; %m^2/s

%% By-pass the preceeding and load existing file
load apef_thorpe_wholeTM.mat


%% Reorder vertical profiles (compute N)
disp('Reordering profiles')
windowRho = Ritp1;
windowTime = timeSensor1;
H = size(windowRho,1);

% N2 calculation
[windowRhoSort, I] = sort(windowRho,1);
[dRx, dRz] = gradient(windowRhoSort);
N2bg = g./windowRhoSort.*dRz./dz;

% Compute local N2
dMat = nan(length(zVecReg), length(windowTime));
N2Mat = nan(length(zVecReg), length(windowTime));
RaVec = []; RaVec2 = [];
for j = 1:length(windowTime)
    rhoVec_raw =  windowRho(:,j);
    [rhoVec_sort, I] = sort(rhoVec_raw);
    d = zVecReg - zVecReg(I);
    rho0 = nanmean(windowRho(:,j));
    RaVec = [RaVec2 g./H/rho0.*(EP1-EP2).*H.^3./nu_cst./kappa_T];
    RaVec2 = [RaVec2 g./H/rho0.*(EP1-EP2).*rms(d).^3./nu_cst./kappa_T];
    dMat(:,j) = d';
    
    % local stratification
    N2local = g./rho0.*gradient(rhoVec_sort)./gradient(zVecReg');
    N2Mat(:,j) = N2local;
end


%% Vertically-averaged fields
time2 = time2maxV('./roc12.mat',timeSensor); 
GamVec = JbVec./epsVec;
N2Vec = nanmean(N2Mat,1);
KVec = GamVec.*epsVec./N2Vec;
RaVec = RaVec2; % <------------- Chose which Ra here!!!
TNVec = (maxdVec.^2./JbVec).^(1/3);
S2bgVec = nanmean(S2bg,1);

%% Average temperature field relative to max vel
dtAve = .05;
timeAve = [-.5:dtAve:.5];
epsAve = nan(1,length(timeAve)); eps2p5 = nan(1,length(timeAve)); eps95p5 = nan(1,length(timeAve));
JbAve = nan(1,length(timeAve)); Jb2p5 = nan(1,length(timeAve)); Jb97p5 = nan(1,length(timeAve));
N2Ave = nan(1,length(timeAve)); N22p5 = nan(1,length(timeAve)); N297p5 = nan(1,length(timeAve));
KAve = nan(1,length(timeAve)); K2p5 = nan(1,length(timeAve)); K97p5 = nan(1,length(timeAve));
RaAve = nan(1,length(timeAve)); Ra2p5 = nan(1,length(timeAve)); Ra97p5 = nan(1,length(timeAve));
TNAve = nan(1,length(timeAve)); TN2p5 = nan(1,length(timeAve)); TN97p5 = nan(1,length(timeAve));
%N2bgAve = nan(1,length(timeAve)); N2bg2p5 = nan(1,length(timeAve)); N2bg97p5 = nan(1,length(timeAve));
S2bgAve = nan(1,length(timeAve)); S2bg2p5 = nan(1,length(timeAve)); S2bg97p5 = nan(1,length(timeAve));

for i = 1:length(timeAve)
   I = find(time2>=timeAve(i)-dtAve & time2<=timeAve(i)+dtAve);
   epsAve(i) = nanmean(epsVec(I));
   JbAve(i) = nanmean(JbVec(I));
   GamAve(i) = nanmean(GamVec(I));   
   N2Ave(i) = nanmean(N2Vec(I));
   KAve(i) = nanmean(KVec(I));
   RaAve(i) = nanmean(RaVec(I));
   TNAve(i) = nanmean(TNVec(I));
   %N2bgAve(i) = nanmean(N2bgVec(I));
   S2bgAve(i) = nanmean(S2bgVec(I));

   %% bootstrap
   nboot = 500;
   N = length(I);
   for iboot = 1:nboot
       % Log-Normal distrib. bootstrap
       r = rand(N,1);
       r = ceil(r*N/1);
       
       eps_boot_b(iboot) = nanmean(epsVec(I(r)));
       Jb_boot_b(iboot) = nanmean(JbVec(I(r)));
       Gam_boot_b(iboot) = nanmean(GamVec(I(r)));
       N2_boot_b(iboot) = nanmean(N2Vec(I(r)));
       K_boot_b(iboot) = nanmean(KVec(I(r)));
       Ra_boot_b(iboot) = nanmean(RaVec(I(r)));
       TN_boot_b(iboot) = nanmean(TNVec(I(r)));
       %      N2bg_boot_b(iboot) = nanmean(N2bgVec(I(r)));
       S2bg_boot_b(iboot) = nanmean(S2bgVec(I(r)));

   end
 
   CI_2p5 = round(2.5/100*nboot);
   CI_97p5 = round(97.5/100*nboot);

   epsSort = sort(eps_boot_b);
   JbSort = sort(Jb_boot_b);
   GamSort = sort(Gam_boot_b);
   N2Sort = sort(N2_boot_b);
   KSort = sort(K_boot_b);
   RaSort = sort(Ra_boot_b);
   TNSort = sort(TN_boot_b);
   %   N2bgSort = sort(N2bg_boot_b);
   S2bgSort = sort(S2bg_boot_b);

   eps2p5(i) = epsSort(CI_2p5);
   eps97p5(i) = epsSort(CI_97p5); 
   Jb2p5(i) = JbSort(CI_2p5);
   Jb97p5(i) = JbSort(CI_97p5); 
   Gam2p5(i) = GamSort(CI_2p5);
   Gam97p5(i) = GamSort(CI_97p5);    
   N22p5(i) = N2Sort(CI_2p5);
   N297p5(i) = N2Sort(CI_97p5); 
   K2p5(i) = KSort(CI_2p5);
   K97p5(i) = KSort(CI_97p5); 
   Ra2p5(i) = RaSort(CI_2p5);
   Ra97p5(i) = RaSort(CI_97p5); 
   TN2p5(i) = TNSort(CI_2p5);
   TN97p5(i) = TNSort(CI_97p5);
   %   N2bg2p5(i) = N2bgSort(CI_2p5);
   %   N2bg97p5(i) = N2bgSort(CI_97p5);
   S2bg2p5(i) = S2bgSort(CI_2p5);
   S2bg97p5(i) = S2bgSort(CI_97p5);
end   

%% Few plots
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 4; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.12; % very left of figure
rigs = 0.14; % very right of figure
tops = 0.03; % top of figure
bots = 0.08; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

x = [-.1 .1 .1 -.1 -.1];
x2 = timeAve;

figure(1)
clf
set(gcf, 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 16 16])

%% S1
s1 = subplot(4,1,1);
YLIM = [1e-8 1e-6];
patch(x, [YLIM(1) YLIM(1) YLIM(2) YLIM(2) YLIM(1)], [1 1 1]*.8, 'linestyle', 'none');
set(gca, 'yscale', 'log')
hold on
%YLIM = get(gca, 'ylim')
%patch(x, [YLIM(1) YLIM(1) YLIM(2) YLIM(2) YLIM(1)], [1 1 1]*.8, 'linestyle', 'none');

yy1 = eps2p5;
yy2 = eps97p5;
patch([x2 fliplr(x2) x2(1)], [yy1 fliplr(yy2) yy1(1)], [.5 .5 .5], 'linestyle', 'none')
semilogy(timeAve, epsAve, 'color', [.5 .5 .5], 'linewidth', 2)
semilogy(timeAve, eps2p5, 'color', [.5 .5 .5], 'linewidth', 2)
semilogy(timeAve, eps97p5, 'color', [.5 .5 .5], 'linewidth', 2)

yy1 = Jb2p5;
yy2 = Jb97p5;
patch([x2 fliplr(x2) x2(1)], [yy1 fliplr(yy2) yy1(1)], [0 0 0], 'linestyle', 'none')
semilogy(timeAve, JbAve, 'color', [.5 .5 .5]*0, 'linewidth', 2)
semilogy(timeAve, Jb2p5, 'color', [.5 .5 .5]*0, 'linewidth', 2)
semilogy(timeAve, Jb97p5, 'color', [.5 .5 .5]*0, 'linewidth', 2)
ylim(YLIM);
set(gca, 'xticklabel', '')
set(gca, 'box', 'on')
set(gca, 'tickdir', 'out')
set(gca, 'xtick', [-.5:.2:.5])
set(gca, 'xgrid', 'on')
ylabel('\epsilon, J_b (W kg^{-1})', 'fontSize', FS)
hold off
text(.45, 1.5e-8, 'a', 'fontSize', FS2, 'fontWeight', 'bold')
text(.52, 1e-7, 'J_b', 'fontSize', FS2, 'fontWeight', 'bold')
text(.52, 8e-7, '\epsilon', 'fontSize', FS2, 'fontWeight', 'bold')
set(gca, 'fontSize', FS)
adjust_space
box on


%% S2
s2 = subplot(4,1,2);
YLIM = [0 .5];
patch(x, [YLIM(1) YLIM(1) YLIM(2) YLIM(2) YLIM(1)], [1 1 1]*.8, 'linestyle', 'none');
hold on
yy1 = Gam2p5;
yy2 = Gam97p5;
patch([x2 fliplr(x2) x2(1)], [yy1 fliplr(yy2) yy1(1)], [0 0 0], 'linestyle', 'none')
ylim(YLIM);
plot(timeAve, GamAve, 'k', 'linewidth', 2)
plot(timeAve, Gam2p5, 'k', 'linewidth', 2)
plot(timeAve, Gam97p5, 'k', 'linewidth', 2)
%hold on
plot(timeAve, GamAve*0+.2, '--k', 'linewidth', 1)
ylabel('\Gamma', 'fontSize', FS)
set(gca, 'xticklabel', '')
set(gca, 'box', 'on')
set(gca, 'tickdir', 'out')
set(gca, 'xtick', [-.5:.2:.5])
set(gca, 'xgrid', 'on')
text(.45, 0.05, 'b', 'fontSize', FS2, 'fontWeight', 'bold')
set(gca, 'fontSize', FS)
adjust_space
box on

%% S3
s3 = subplot(4,1,3);
YLIM = [2e10 5e12];
patch(x, [YLIM(1) YLIM(1) YLIM(2) YLIM(2) YLIM(1)], [1 1 1]*.8, 'linestyle', 'none');
set(gca, 'yscale', 'log')

hold on
yy1 = Ra2p5;
yy2 = Ra97p5;
patch([x2 fliplr(x2) x2(1)], [yy1 fliplr(yy2) yy1(1)], [0 0 0], 'linestyle', 'none')
ylim(YLIM);
set(gca, 'xticklabel', '')
set(gca, 'box', 'on')
set(gca, 'tickdir', 'out')
set(gca, 'xtick', [-.5:.2:.5])
set(gca, 'xgrid', 'on')
ylabel('Ra', 'fontSize', FS)
semilogy(timeAve, RaAve, 'k', 'linewidth', 2)
semilogy(timeAve, Ra2p5, 'k', 'linewidth', 2)
semilogy(timeAve, Ra97p5, 'k', 'linewidth', 2)
hold off
text(.45, 3e10, 'c', 'fontSize', FS2, 'fontWeight', 'bold')
set(gca, 'fontSize', FS)
adjust_space
box on

% $$$ %% S3
% $$$ s3 = subplot(4,1,3);
% $$$ YLIM = [1e-6 1e-4];
% $$$ patch(x, [YLIM(1) YLIM(1) YLIM(2) YLIM(2) YLIM(1)], [1 1 1]*.8, 'linestyle', 'none');
% $$$ set(gca, 'yscale', 'log')
% $$$ 
% $$$ hold on
% $$$ yy1 = N22p5;
% $$$ yy2 = N297p5;
% $$$ patch([x2 fliplr(x2) x2(1)], [yy1 fliplr(yy2) yy1(1)], [0 0 0], 'linestyle', 'none')
% $$$ ylim(YLIM);
% $$$ set(gca, 'xticklabel', '')
% $$$ set(gca, 'box', 'on')
% $$$ set(gca, 'tickdir', 'out')
% $$$ set(gca, 'xtick', [-.5:.2:.5])
% $$$ set(gca, 'xgrid', 'on')
% $$$ ylabel('N^2 (s^{-2})', 'fontSize', FS)
% $$$ semilogy(timeAve, N2Ave, 'k', 'linewidth', 2)
% $$$ semilogy(timeAve, N22p5, 'k', 'linewidth', 2)
% $$$ semilogy(timeAve, N297p5, 'k', 'linewidth', 2)
% $$$ hold off
% $$$ text(.45, 1.5e-6, 'c', 'fontSize', FS2, 'fontWeight', 'bold')
% $$$ set(gca, 'fontSize', FS)
% $$$ adjust_space
% $$$ box on

%% S4
s4 = subplot(4,1,4);
YLIM = [1e-3 1e-1];
patch(x, [YLIM(1) YLIM(1) YLIM(2) YLIM(2) YLIM(1)], [1 1 1]*.8, 'linestyle', 'none');
set(gca, 'yscale', 'log')
set(gca, 'box', 'on')
set(gca, 'tickdir', 'out')
set(gca, 'xtick', [-.5:.2:.5])
set(gca, 'xgrid', 'on')
hold on
yy1 = K2p5;
yy2 = K97p5;
patch([x2 fliplr(x2) x2(1)], [yy1 fliplr(yy2) yy1(1)], [0 0 0], 'linestyle', 'none')
ylim(YLIM);
semilogy(timeAve, KAve, 'k', 'linewidth', 2)
semilogy(timeAve, K2p5, 'k', 'linewidth', 2)
semilogy(timeAve, K97p5, 'k', 'linewidth', 2)
ylabel('K (m^2 s^{-1})', 'fontSize', FS)
xlabel('\phi (day)', 'fontSize', FS)
text(.45, 1.5e-3, 'd', 'fontSize', FS2, 'fontWeight', 'bold')
set(gca, 'fontSize', FS)
adjust_space
box on


I = find(timeAve>=-.1 & timeAve<= .1);
II = find(timeAve<-.1 | timeAve> .1);

[nanmean(epsAve(I)) nanmean(eps2p5(I)) nanmean(eps97p5(I))]
[nanmean(epsAve(II)) nanmean(eps2p5(II)) nanmean(eps97p5(II))]
[nanmean(JbAve(I)) nanmean(Jb2p5(I)) nanmean(Jb97p5(I))]
[nanmean(JbAve(II)) nanmean(Jb2p5(II)) nanmean(Jb97p5(II))]
[nanmean(KAve(I)) nanmean(K2p5(I)) nanmean(K97p5(I))]
[nanmean(KAve(II)) nanmean(K2p5(II)) nanmean(K97p5(II))]
[nanmean(N2Ave(I)) nanmean(N22p5(I)) nanmean(N297p5(I))]
[nanmean(N2Ave(II)) nanmean(N22p5(II)) nanmean(N297p5(II))]
[nanmean(GamAve(I)) nanmean(Gam2p5(I)) nanmean(Gam97p5(I))]
[nanmean(GamAve(II)) nanmean(Gam2p5(II)) nanmean(Gam97p5(II))]

set(gcf, 'renderer', 'painters')
outfile = 'RaTidal.eps';
print('-depsc2', outfile)


%% Manuscript revision plot (second version of the figure with convection time scale)
epsK_plots_plot


epsK_plots2