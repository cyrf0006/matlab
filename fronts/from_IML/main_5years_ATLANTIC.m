%   Guillaume Desbiens le 28 novembre 2011
%
%   Calcul des moyennes des probabilit�s de presence de fronts
%   Ce programme va chercher les matrices de probabilit�s mensuelles de d�tection de 
%   fronts et va calculer la climatologie entre 2006 et 2010
%   Les matrices sources se trouvent sous:J:\matlab_atlantic\moyenne_INPUT\
%   Le r�sulat sera enregistr� sous J:\matlab_atlantic\probabilite_OUTPUT\annuelle
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
    source ='J:\matlab_atlantic\moyenne_INPUT'; 
%------------------------------------------------------------------------

   % le dossier de destination DOIT �TRE CR�� AVANT l'�x�cution 

    dest ='J:\matlab_atlantic\probabilite_OUTPUT\annuelle\';                                   
   
[mean_prob_5years,cumul_prob,used_files60]=prob_5years_ATLANTIC();
prob_5years_integer=uint8(mean_prob_5years); % transforme la matrice mean_prob_25years en integer (4.4 devient 4)
                                                % n�cessaire pour le fichier HDF
    %------------------------
    % FIGURE
    %--------------------------

    h=figure(1);
    imagesc(lon(1,:)',lat(:,1),mean_prob_5years)

    set(gca,'YDir','normal')

    patch(coast(1:end-1,1),coast(1:end-1,2),[0.75 0.75 0.75]) % ajoute le trait de c�te � l'image

    caxis([0 20])   % Echelle de couleur 0 a 20% de probabilit� (les valeurs >20 sont toutes de la m�me couleur)
    
    colormap(mymap)
    colorbar


    title(' Front detection probability averaged over 2006-2010 (%)') 
                                                                           
    save([dest 'Prob_2006-2010.mat'],'lat','lon','used_files60','mean_prob_5years','cumul_prob')
    saveas(h,[dest 'Prob_2006-2010'],'png') 
    saveas(h,[dest 'Prob_2006-2010'],'fig')

    cd('J:\matlab_atlantic\probabilite_OUTPUT\HDF\yearly\');  
    imwrite(prob_5years_integer,'Prob_2006-2010.hdf');  % cr�e le fichier HDF
    
    delete(h)

%  end       