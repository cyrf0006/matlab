%   Guillaume Desbiens le 15 d�cembre 2011
%
%   Calcul l'�cart-type des probabilit�s de presence de fronts
%   Ce programme va chercher les matrices de probabilit�s mensuelles de d�tection de 
%   fronts et va calculer l'�cart-type pour chaque pixel pour la p�riode 1986 et 2010
%   Les matrices sources se trouvent sous:J:\matlab_atlantic\moyenne_INPUT\
%   Le r�sulat sera enregistr� sous
%   J:\matlab_atlantic\probabilite_OUTPUT\stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

clc;
clear all;
close all;

dim_im=[3499 3000]; % taille de l'image 

    
%Pour mettre le fond blanc dans l'image
mymap=colormap(jet(25));
mymap(1,:)=[0.9 0.9 0.9];

%-----Caract�ristiques de la r�gion atlantic---------------
%-----LATITUDE/LONGITUDE------------------------------%

info = hdfinfo('atlantic_latlon.hdf');
dsets = info.SDS;
lat=hdfread(dsets(1));
lon=hdfread(dsets(2));


%-----TRAIT DE COTE------------------------------%

coast=load('coastline_atlantic.dat');
infoland=hdfinfo('landmask_atlantic_binary.hdf');
dsetland=infoland.SDS;
land = hdfread(dsetland);
land=double(land);

%------------------------------------------------------------------------
    source ='E:\matlab_atlantic\moyenne_INPUT'; 
%------------------------------------------------------------------------

   % le dossier de destination DOIT �TRE CR�� AVANT l'�x�cution 

    dest ='E:\matlab_atlantic\probabilite_OUTPUT\stats\';                                   
   
[mean_gamma_25years,cumul_gamma,used_files300]=stats_25years_ATLANTIC();
% prob_25years_integer=uint8(mean_prob_25years); % transforme la matrice mean_prob_25years en integer (4.4 devient 4)
                                                % n�cessaire pour le fichier HDF
                                               
                                                
    %------------------------
    % FIGURE
    %--------------------------

    h=figure(1);
    imagesc(lon(1,:)',lat(:,1),mean_prob_25years)

    set(gca,'YDir','normal')

    patch(coast(1:end-1,1),coast(1:end-1,2),[0.75 0.75 0.75]) % ajoute le trait de c�te � l'image

    caxis([0 20])   % Echelle de couleur 0 a 20% de probabilit� (les valeurs >20 sont toutes de la m�me couleur)

    
    colormap(mymap)
    colorbar


    title(' Front detection probability averaged over 1986-2010 (%)') 
    
                                                                        
    save([dest 'Prob_1986-2010.mat'],'lat','lon','used_files300','mean_prob_25years','cumul_prob')
    saveas(h,[dest 'Prob_1986-2010'],'png') 
    saveas(h,[dest 'Prob_1986-2010'],'fig')

    cd('E:\matlab_atlantic\probabilite_OUTPUT\HDF\yearly\');  
    imwrite(prob_25years_integer,'Prob_1986-2010.hdf');  % cr�e le fichier HDF
    
    delete(h)
