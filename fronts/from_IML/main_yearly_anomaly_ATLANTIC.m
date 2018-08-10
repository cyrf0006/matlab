%   Guillaume Desbiens le 15 d�cembre 2011
%
%   Cr�ation des images annuelles d'anomalie de probabilit�s de d�tection de fronts
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
    source ='J:\matlab_atlantic\moyenne_INPUT\'; 
%------------------------------------------------------------------------

% le dossier de destination DOIT �TRE CR�� AVANT l'�x�cution 
dest ='J:\matlab_atlantic\probabilite_OUTPUT\anomaly\';  
    
extension='.mat';

%Ouverture du fichier moyenne annuelle totale 
data=load('J:\matlab_atlantic\probabilite_OUTPUT\annuelle\Prob_1986-2010.mat');

%lecture de la probabilit� moyenne totale de toutes les ann�es
mean_25y=data.mean_prob_25years;

%---calcul et cr�ation d'une image d'anomalie par ann�e
for y=2000
        year=num2str(y);
%     [STD_sigma_year,sigma_year,sigma_25years,cumul_sigma,used_files25]=anomaly_year_ATLANTIC(source,dim_im,extension);
    [STD_sigma_year,sigma_year,used_files25]=anomaly_year_ATLANTIC();
    
                                               
                                                
    %------------------------
    % FIGURE
    %--------------------------

    h=figure(1);
    imagesc(lon(1,:)',lat(:,1),STD_sigma_year)

    set(gca,'YDir','normal')

    patch(coast(1:end-1,1),coast(1:end-1,2),[0.75 0.75 0.75]) % ajoute le trait de c�te � l'image

    caxis([0 10])   % Echelle de couleur 0 a 20% de probabilit� (les valeurs >20 sont toutes de la m�me couleur)

    
    colormap(mymap)
    colorbar


    title([' Standardized anomaly of front detection for' year]) 
    
                                                                        
    save([dest 'Anomaly_' year '.mat'],'STD_sigma_year','sigma_year')
    saveas(h,[dest 'Anomaly_' year ],'png') 
    saveas(h,[dest 'Anomaly_' year ],'fig')

%     cd('J:\matlab_atlantic\probabilite_OUTPUT\HDF\anomaly\');  
%     imwrite(prob_25years_integer,'Prob_1986-2010.hdf');  % cr�e le fichier HDF
%     clear ('STD_sigma_year','sigma_year'); 
    delete(h)
end
