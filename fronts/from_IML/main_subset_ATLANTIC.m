% Ce script va chercher les images recherchées selon l'intervalle de temps
% désiré (J:\matlab_atlantic\probabilite_OUTPUT\mensuelle) et créer UNE 
% image de probabilités de présence de fronts  PAR MOIS dans le répertoire ('probabilite_OUTPUT\mensuel\').

clc
clear all
close all
dim_im=[3499 3000]; % taille de l'image 

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

%--------------------------------------------------------------------------
% Spécifier l'intervalle temporel de recherche des images
%--------------------------------------------------------------------------
 
for y=2000  % A MODIFIER
    year=num2str(y);
    for m=1    % A MODIFIER selon les mois à traiter
        month=m;
        if (month<10)        
            strmonth= ['0' num2str(month)];
        else
            strmonth= num2str(month);
        end
        
       %-------------------------------------------------------------------
        source = 'J:\matlab_atlantic\moyenne_INPUT'; 
       %------------------------
       
       
       % le répertoire de destination DOIT ÊTRE CRÉÉ AVANT l'éxécution 
       dest ='J:\matlab_atlantic\subset_output\';                                   
       
        [count_edge_subset]=subset_ATLANTIC(source,year,strmonth,dim_im); 
        
        %count_edge(count_edge ==0)=land(count_edge==0)+count_edge(count_edge==0);  %#ok<AGROW> % pour mettre trait de cote
        
        h=figure(1);
        %imagesc(lon(1,:)',lat(:,1),count_edge)
        imagesc(lon(1,:)',lat(:,1),probability)
        set(gca,'YDir','normal')
        
        patch(coast(1:end-1,1),coast(1:end-1,2),[0.75 0.75 0.75]) % ajoute le trait de côte à l'image

        caxis([0 20])   % Echelle de couleur 0 a 20% de probabilité (les valeurs >20 sont toutes de la même couleur)
        colormap(mymap)
        colorbar('location','southoutside')
        
        
        title(['Monthly front detection probability (%)  ' dest(30:end)]) 
%         title(['Prob. mensuelle détection de front (%) ' dest(30:end) '(' num2str(c) ' jours)']) 
        % dest= ex: probabilite_OUTPUT\mensuel_1995_04 => dest(30:end) conserve donc SEUL. yyyy_mm 
        % "c" provient de GetPict, c'est le nombre de fichiers (jours) utiliser pour le calcul 

        
        save([dest 'Prob_Mensuel_' year '_' strmonth '.mat'],'mymap','cloud_prob','image_overlay','count_edge','count_pix','probability')                                                                            
%         save([dest2 'Prob_Mensuel_' year '_' strmonth '.mat'],'cloud_prob','image_overlay','count_edge','count_pix','c','probability')

        
        saveas(h,[dest 'Prob_'  year '_' strmonth],'png') 
        saveas(h,[dest 'Prob_' year '_' strmonth],'fig')
 
        
        delete(h)
    end
end
