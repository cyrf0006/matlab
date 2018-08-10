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
load grid.mat;

timeVec = time(:,3);
[Y,I] = min(abs(timeVec-datenum(2000,1,6,6,0,0)));


% $$$ mooringDepth = 919;
% $$$ noSteps = size(time,1);
% $$$ I = find(zeta>=794 & zeta <= 919); % mooring depth range
% $$$ Iu = find(z>=794 & z <= 919); % mooring depth range
% $$$ [Y, xI] = min(abs(H0-mooringDepth));
% $$$ XLIM = [datenum(2000,1,5,19.5,0,0) datenum(2000,1,6,19.5,0,0)];
% $$$ 
% $$$ 
% $$$ zVec = zeta(I);
% $$$ zuVec = z(Iu);
% $$$ sigMat = nan(length(I), noSteps);
% $$$ uMat = nan(length(Iu), noSteps);
% $$$ vMat = nan(length(Iu), noSteps);
% $$$ timeVec = time(:,3);


sig = getfield('sig',I);
sig = putNaN(sig,flag_s);
sig = sig';

figure(1)
contourf(x, zeta, sig, 100, 'lineStyle', 'none')      
hold on
contour(x, zeta, sig, 50, 'color', 'k')    


% $$$ 
% $$$ u = getfield('u',i);
% $$$ u = putNaN(u,flag_u);
% $$$ u = u';    
% $$$ v = getfield('v',i);
% $$$ v = putNaN(v,flag_u);
% $$$ v = v';      

