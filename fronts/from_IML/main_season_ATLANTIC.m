% Ce script va chercher les images recherchées selon l'intervalle de temps
% désiré (J:\matlab_atlantic\moyenne_INPUT\yyyy) et créer UNE CARTE de 
% probabilités de présence de fronts  PAR SAISON dans les répertoires:
%  (J://probabilite_OUTPUT\saisonniere\')

% IMPORTANT: si les saisons sont divisées telles que par Environnement Canada
% les lignes 42 à 51 et 84 à 92 doivent être mises en commentaires alors que
% 55 à 64 et 95 à 103 doivent être décommentées,

clc
clear all
close all
dim_im=[3499 3000]; % taille de l'image 

season=4;           %A MODIFIER: 1 pour hiver, 2 printemps, 3 été, 4 automne

%Pour mettre le fond blanc dans l'image
mymap=colormap(jet(25));
mymap(1,:)=[0.9 0.9 0.9];
% mymap(2:4,:)=[1 1 1;1 1 1;1 1 1];


%-----Caractéristiques de la région gulf---------------
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
%land(land==1)=-1;   % la terre est déja à -1 (255 dans le hdf)À VÉRIFIER
land=double(land);

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


% DÉCOMMENTER pour SAISONS SELON ENVIRONNEMENT CANADA
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


%--------------------------------------------------------------------------
% Spécifier l'intervalle temporel de recherche des images
%--------------------------------------------------------------------------

for y=1986:2010  % A MODIFIER
    year=num2str(y);
   
    
    source ='J:\matlab_atlantic\moyenne_INPUT'; 

   %------------------------

   % le dossier de destination DOIT ÊTRE CRÉÉ AVANT l'éxécution 

    dest ='J:\matlab_atlantic\probabilite_OUTPUT\saisonniere\';                                   
   
   % COMMENTER pour SAISONS SELON ENVIRONNEMENT CANADA
%     if (season==1)        
%            [prob_season,cumul_prob,used_files3]=prob_winter_ATLANTIC(year);
%     elseif (season==2)
%            [prob_season,cumul_prob,used_files3]=prob_spring_ATLANTIC(year);
%     elseif (season==3)
%            [prob_season,cumul_prob,used_files3]=prob_summer_ATLANTIC(year);
%     elseif (season==4)
%            [prob_season,cumul_prob,used_files3]=prob_autumn_ATLANTIC(year);
%     end
    
    % DÉCOMMENTER pour SAISONS SELON ENVIRONNEMENT CANADA
    if (season==1)        
           [prob_season,cumul_prob,used_files3]=ECprob_winter_ATLANTIC(year);
    elseif (season==2)
           [prob_season,cumul_prob,used_files3]=ECprob_spring_ATLANTIC(year);
    elseif (season==3)
           [prob_season,cumul_prob,used_files3]=ECprob_summer_ATLANTIC(year);
    elseif (season==4)
           [prob_season,cumul_prob,used_files3]=ECprob_autumn_ATLANTIC(year);
    end
    
    %------------------------
    % FIGURE
    %--------------------------

    h=figure(1);
    imagesc(lon(1,:)',lat(:,1),prob_season)

    set(gca,'YDir','normal')

    patch(coast(1:end-1,1),coast(1:end-1,2),[0.75 0.75 0.75]) % ajoute le trait de côte à l'image

    caxis([0 20])   % Echelle de couleur 0 a 20% de probabilité (les valeurs >20 sont toutes de la même couleur)
    
    colormap(mymap)
    colorbar


    title([ strseason '_ ' year ' front detection probability (%)'  ]) 
                                                                               
    save([dest 'Prob_' strseason '_' year '.mat'],'lat','lon','used_files3','prob_season','cumul_prob')
    saveas(h,[dest 'Prob_' strseason '_' year ],'png') 
    saveas(h,[dest 'Prob_' strseason '_' year ],'fig')
  
    delete(h)
    clear('prob_season','cumul_prob') % permet à MATLAB de traiter plusieurs années en boucle sans manquer de mémoire
                                      % COMMENTER pour que les variables s'affichent 
                                      % (valider au moins un fichier, la variable used_files3 devrait être ==3 (3 mois * 1 an) 
end
