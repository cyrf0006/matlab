%   Guillaume Desbiens le 28 novembre 2011
%
%   Calcul des moyennes des probabilit�s de presence de fronts
%   Ce programme va chercher les matrices de probabilit�s mensuelles de d�tection de 
%   fronts et va calculer les climatologies selon la saison choisie 
%   Les matrices sources se trouvent sous:J:\matlab_atlantic\moyenne_INPUT\
%   Les r�sulats sont enregistr�s sous J:\matlab_atlantic\probabilite_OUTPUT\saisonniere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% IMPORTANT: si les saisons sont divis�es telles que par Environnement Canada
% les lignes 22-31 et 79-87 doivent �tre mises en commentaires alors que
% 34-43 et 90-98 doivent �tre d�comment�es,

clc;
clear all;
close all;

season=4;           %A MODIFIER: 1 pour hiver, 2 printemps, 3 �t�, 4 automne

dim_im=[3499 3000]; % taille de l'image 

% COMMENTER pour SAISONS SELON ENVIRONNEMENT CANADA
% if (season==1)        
%         strseason='winter';
%     elseif (season==2)
%         strseason='spring';
%     elseif (season==3)
%         strseason='summer';
%     elseif (season==4)
%         strseason='autumn';
%     else disp ('ERREUR: la saison choisie n est pas valide');    
% end

% D�COMMENTER pour SAISONS SELON ENVIRONNEMENT CANADA
if (season==1)        
        strseason='winter(Dec-Feb)';
    elseif (season==2)
        strseason='spring(March-May)';
    elseif (season==3)
        strseason='summer(June-Aug)';
    elseif (season==4)
        strseason='autumn(Sept-Nov)';
    else disp ('ERREUR: la saison choisie n est pas valide');    
end

    
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

%------------------------------------------------------------------------
source ='J:\matlab_atlantic\moyenne_INPUT'; 

   %------------------------

   % le dossier de destination DOIT �TRE CR�� AVANT l'�x�cution 

dest ='J:\matlab_atlantic\probabilite_OUTPUT\saisonniere\';                                   
   

% COMMENTER pour SAISONS SELON ENVIRONNEMENT CANADA
% if (season==1)        
%        [mean_prob_season,cumul_prob,used_files15]=prob_winter_5years_ATLANTIC();
% elseif (season==2)
%        [mean_prob_season,cumul_prob,used_files15]=prob_spring_5years_ATLANTIC();
% elseif (season==3)
%        [mean_prob_season,cumul_prob,used_files15]=prob_summer_5years_ATLANTIC();
% elseif (season==4)
%        [mean_prob_season,cumul_prob,used_files15]=prob_autumn_5years_ATLANTIC();
% end

% D�COMMENTER pour SAISONS SELON ENVIRONNEMENT CANADA
if (season==1)        
       [mean_prob_season,cumul_prob,used_files15]=ECprob_winter_5years_ATLANTIC();
elseif (season==2)
       [mean_prob_season,cumul_prob,used_files15]=ECprob_spring_5years_ATLANTIC();
elseif (season==3)
       [mean_prob_season,cumul_prob,used_files15]=ECprob_summer_5years_ATLANTIC();
elseif (season==4)
       [mean_prob_season,cumul_prob,used_files15]=ECprob_autumn_5years_ATLANTIC();
end

    %------------------------
    % FIGURE
    %--------------------------

    h=figure(1);
    imagesc(lon(1,:)',lat(:,1),mean_prob_season)

    set(gca,'YDir','normal')

    patch(coast(1:end-1,1),coast(1:end-1,2),[0.75 0.75 0.75]) % ajoute le trait de c�te � l'image

    caxis([0 20])   % Echelle de couleur 0 a 20% de probabilit� (les valeurs >20 sont toutes de la m�me couleur)
    
    colormap(mymap)
    colorbar


    title([ strseason ' front detection probability (%) averaged over 2006-2010'  ]) 
                                                                      
    save([dest 'Prob_' strseason '_2006-2010.mat'],'lat','lon','used_files15','mean_prob_season','cumul_prob')
    saveas(h,[dest 'Prob_'  strseason '_2006-2010'],'png') 
    saveas(h,[dest 'Prob_' strseason '_2006-2010'],'fig')
  
    delete(h)
      