clear all
close all

% $$$ VIZ: La commande pour visualizer rapidement un champ. Exemple: 
% $$$ 
% $$$ viz('u',12);
% $$$ 
% $$$ où 12 est le frame de output.  viz('u',1) donne la solution initiale. 
% $$$ 
% $$$ GETFIELD: La commande pour extraire les données (elle est utilisée dans VIZ). Exemple:
% $$$ 
% $$$ sig = getfield('sig',1100);
% $$$ h = getfield('h',1100);  % Pour l'élévation de la surface
% $$$ 
% $$$ PUTNAN: Pour mettre des NaN sous la topographie. Exemple
% $$$ 
% $$$ sig = putNaN(sig,flag_s);
% $$$ u    = putNaN(u,flag_u);
% $$$ v    = putNaN(v,flag_u);
% $$$ w    = putNaN(w,flag_s);
% $$$ 
% $$$ La grille est déjà en Matlab dans le fichier grid.mat
% $$$ x: Les points de grille de u et sigma
% $$$ chi: les point de grille de w
% $$$ z: les u
% $$$ zeta: les sigma et w
% $$$ 
% $$$ flag_s: le flag des scalaires (mais aussi de w)
% $$$ flag_u: le flag de u et v

%% some info on the grid
load time.out;
load grid.mat
mooringDepth = 919;

timeVec = time(:,3);
XLIM = [datenum(2000,1,5,19.5,0,0) datenum(2000,1,6,19.5,0,0)];
XLIM = [datenum(2000,1,6,0,0,0) datenum(2000,1,6,8,0,0)];

Itime = find(timeVec>=XLIM(1) & timeVec<= XLIM(end));
noSteps = length(Itime);
timeVec = timeVec(Itime);

I = find(zeta>=794 & zeta <= 919); % mooring depth range
Iu = find(z>=794 & z <= 919); % mooring depth range
[Y, xI] = min(abs(H0-mooringDepth));

zVec = zeta(I);
zuVec = z(Iu);
sigMat = nan(length(I), noSteps);
uMat = nan(length(Iu), noSteps);
vMat = nan(length(Iu), noSteps);
NMat = nan(length(Iu), noSteps);


disp('Extract model results...')
for i = 1:noSteps
    if mod(i, 100) == 0;
        disp(sprintf('   timestep %d/%d', i, noSteps))
    end
    sig = getfield('sig',Itime(i));
    sig = putNaN(sig,flag_s);
    sig = sig';
    u = getfield('u',Itime(i));
    u = putNaN(u,flag_u);
    u = u';    
    v = getfield('v',Itime(i));
    v = putNaN(v,flag_u);
    v = v';      
    
    sigMat(:,i) = sig(I, xI);        
    uMat(:,i) = u(Iu, xI);        
    vMat(:,i) = v(Iu, xI);        
    rhoVec = sort(sig(I, xI));
    NVec = sqrt(9.8./1025.*gradient(rhoVec,1,1)./gradient(zVec',1,1));
    NMat(:,i) = NVec;
end
disp('done!')



%%%%%%%% ------------ Model ---------- %%%%%%%%

% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 2; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.03; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.1; % very right of figure
tops = 0.03; % top of figure
bots = 0.07; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

figure(1)
clf
set(gcf, 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 16 14])
subplot(2,1,1)
contourf(timeVec, zVec, sigMat, 100, 'lineStyle', 'none')
hold on
contour(timeVec, zVec, sigMat, 50, 'color', 'k')
datetick('x',15)
xlim(XLIM)
set(gca, 'xticklabel', []);
c = colorbar;
caxis([31.9 32.1])
adjust_space
pause(.5)
c_pos = get(c, 'pos');
c_pos(3) = c_pos(3)/2;
c_pos(1) = c_pos(1)-0.02;
c_pos(2) = c_pos(2)+.02;
c_pos(4) = c_pos(4)-.04;
set(c, 'pos', c_pos);
set(gca, 'ydir', 'reverse')
text(XLIM(1)+1/100, 900, '\sigma_1(kg m^{-3})', 'fontWeight', 'bold', ...
     'backgroundcolor', [1 1 1])
ylabel('Depth (m)')
xlim(XLIM)

subplot(2,1,2)
contourf(timeVec, zVec, NMat, 100, 'lineStyle', 'none')
% $$$ pcolor(timeVec, zVec, log10(NMat))
% $$$ shading flat
hold on
contour(timeVec, zVec, sigMat, 50, 'color', 'k')
datetick('x',15)
xlim(XLIM)
xlabel('time')
%caxis([-3.2 -3])
caxis([0.8e-3 3.2e-3])
c = colorbar;
set(gca, 'ydir', 'reverse')
adjust_space
pause(.5)
c_pos = get(c, 'pos');
c_pos(3) = c_pos(3)/2;
c_pos(1) = c_pos(1)-0.02;
c_pos(2) = c_pos(2)+.02;
c_pos(4) = c_pos(4)-.04;
xlim(XLIM)


text(XLIM(1)+1/100, 900, 'N (s^{-1})', 'fontWeight', 'bold','backgroundcolor', [1 1 1])

set(c, 'pos', c_pos);

set(gcf, 'renderer', 'painter')
print(gcf, '-depsc', 'aestusNT.eps')