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

% snapshot info 
frameTime = datenum(2000,1,6,7,0,0);
XLIM = [-4e4 2e4];
ZLIM = [0 1600];
XLIM = [-3e4 -1e4];
ZLIM = [250 1000];
XLIM = [-2.85e4 -2.6e4];
ZLIM = [450 550];
% $$$ XLIM = [-2.85e4 -2.6e4];
% $$$ ZLIM = [510 545];
% Get data from this frame
timeVec = time(:,3);
[Y, I] = min(abs(timeVec-frameTime));
sig = getfield('sig',I);
sig = putNaN(sig,flag_s);
sig = sig';
u = getfield('u',I);
u = putNaN(u,flag_u);
u = u';    
v = getfield('v',I);
v = putNaN(v,flag_u);
v = v';   
w = getfield('w',I);
w = putNaN(w,flag_s);
w = w'; 

% reduce the grid
I = find(x>=XLIM(1) & x<=XLIM(2));
J = find(z>=ZLIM(1) & z<=ZLIM(2));
II = find(chi>=XLIM(1) & chi<=XLIM(2));
JJ = find(zeta>=ZLIM(1) & zeta<=ZLIM(2));
xVec = x(I);
zVec = z(J);
chiVec = chi(II);
zetaVec = zeta(JJ);
u = u(J,I);
sig = sig(JJ,I);
w = w(JJ,II);


quiverDec = 3;
figure(1)
clf
set(gcf,'PaperUnits','centimeters', 'paperposition', [0 0 18 14])
pcolor(xVec, zetaVec, sig)
shading flat  
set(gca, 'ydir', 'reverse')
hold on
contour(xVec, zetaVec, sig, 100, 'color', 'k')

Z = zVec(1:quiverDec:end);
X = xVec(1:quiverDec:end);
U = u(1:quiverDec:end, 1:quiverDec:end);
W = w(1:quiverDec:end, 1:quiverDec:end);


[kmax imax] = size(U);
vecColor= [1 1 1];
unit = false;
scale = .25;

for i = 1:imax
    for k = 1:kmax
        arrow7(X(i),Z(k),U(k,i),W(k,i),scale,vecColor,unit);
    end
end

set(gcf, 'renderer', 'painter')
print(gcf, '-depsc', 'aestus_breaking2D.eps')




