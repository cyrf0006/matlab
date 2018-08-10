%   Guillaume Desbiens le 28 novembre 2011
%
%   Calcul des moyennes des probabilit�s de presence de fronts
%   Ce programme va chercher les matrices de probabilit�s mensuelles de d�tection de 
%   fronts et va calculer les climatologies selon l'ann�e choisie 
%   Les matrices sources se trouvent sous:J:\matlab_atlantic\moyenne_INPUT\
%   Les r�sulats sont enregistr�s sous J:\matlab_atlantic\probabilite_OUTPUT\saisonniere
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
%land(land==1)=-1;   % la terre est d�ja � -1 (255 dans le hdf)� V�RIFIER
land=double(land);



for y=2006  % A MODIFIER
    year=num2str(y);

    %------------------------------------------------------------------------
    source ='J:\matlab_atlantic\moyenne_INPUT'; 

       %------------------------

       % le dossier de destination DOIT �TRE CR�� AVANT l'�x�cution 

    dest ='J:\matlab_atlantic\probabilite_OUTPUT\annuelle\'; 

        %     [probability_season,cloud_season,used_files]=probabilities_season_ATLANTIC();
        %     [probability_season,used_files]=probabilities_season_ATLANTIC();
        %     [cloud_season,used_files2]=cloud_season_ATLANTIC();

     [mean_prob_year,cumul_prob,used_files12] = prob_year_ATLANTIC(year);

     prob_year_integer=uint8(mean_prob_year); % transforme la matrice mean_prob_year en integer (4.4 devient 4)
                                                % n�cessaire pour le fichier HDF

            %------------------------
            % FIGURE
            %--------------------------

    h=figure(1);
    imagesc(lon(1,:)',lat(:,1),mean_prob_year)
%     imagesc(lon(1,:)',lat(:,1),cloud_season)

    set(gca,'YDir','normal')

    patch(coast(1:end-1,1),coast(1:end-1,2),[0.75 0.75 0.75]) % ajoute le trait de c�te � l'image

    caxis([0 20])   % Echelle de couleur 0 a 20% de probabilit� (les valeurs >20 sont toutes de la m�me couleur)
%     caxis([0 100])   % Echelle de couleur 0 a 100% de probabilit� 

    colormap(mymap)
    colorbar

    title(['Front detection probability (%) for  ' year ]) 

    save([dest 'Prob_' year '.mat'],'lat','lon','used_files12','mean_prob_year','cumul_prob')
    saveas(h,[dest 'Prob_' year ],'png') 
    saveas(h,[dest 'Prob_' year ],'fig')
    
    
    delete(h)  
    
    clc;
    clear ('cumul_prob','mean_prob_year');
end       