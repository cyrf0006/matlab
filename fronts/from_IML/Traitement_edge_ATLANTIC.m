%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ce script transforme les données de SST en degré Celsius et
% superpose ces images SST coupées et leurs edge . 
% % Entrées :	J:\matlab_atlantic\traitement_INPUT\yyyy
% % •	yyyymmdd...L1..global.msst.mean.1d.2_Cut.hdf
% % •	yyyymmdd...L1..global.msst.mean.1d.2_Cut_Sied.hdf
% % Sorties :	J:\matlab_atlantic\traitement_OUPUT\yyyy
% % •	overlay_yyyymmdd...L1..global.msst.mean.1d.2_Cut..mat  
% % •	overlay_yyyymmdd.fig
% % •	overlay_yyyymmdd.png
%
% La matrice finale .mat contient terre = -5, nuage ou glace = -4, edge  =
% -2.5 et température de -2 à 28 °C.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc 
clear all
close all

%------Pourcentage minimal d'eau sur la carte----
pour_min=2;

mymap=colormap(jet);
mymap(1,:)=[0.9 0.9 0.9];
mymap(2:3,:)=[1 1 1;1 1 1];
mymap(4:5,:)=[1 0 1;1 0 1];


%-----Caractéristiques de la région atlantic---------------
%-----LATITUDE/LONGITUDE------------------------------%

info = hdfinfo('atlantic_latlon.hdf'); % obtenu avec wam_series : prendre n'importe quelle image / save as HDF with LLA
dsets = info.SDS;
lat=double(hdfread(dsets(1))); %s'assurer que le sds 1 soit latitude en l'ouvrant ds wim 
lon=double(hdfread(dsets(2)));
% lat=lat./100;
% lon=lon./100;


%-----TRAIT DE COTE------------------------------%
coast=load('coastline_atlantic.dat'); %obtenu avec le coast extractor via le site de la NOAA  http://www.ngdc.noaa.gov/mgg/coast/
                                      % les ">" doivent être remplacé par "nan nan" 
                                                                                                                                  
infoland=hdfinfo('landmask_atlantic_binary.hdf'); % landmask utilisé pour les traitements sur wam qui a été binariser (terre=255 reste =0)
dsetland=infoland.SDS;
land = hdfread(dsetland);
%land(land==255)=-1; le land est déjà à -1 (Matlab transforme les 255 en -1 quand il load)
land=double(land);



%-----Slope/Intercept---------% proviennent des attributs de n'importe quelle image ds wim
m=0.15; 
b=-3;

for y=2002      % A MODIFIER
    year=num2str(y);
    source=['J:\matlab_atlantic\traitement_INPUT\' year];  % va chercher les données brutes: image SST et image front
    dest=['J:\matlab_atlantic\traitement_OUTPUT\' year '\'];
    
   Dir_source=dir(source);
   dim_source=size(Dir_source);
    
    
    for f=3:2:(dim_source(1)-1)  % 1 et 2 pour les fichiers cachés (.) et (..)
        
        
        %-------------------------%image SST brute%----------------------------------%
        infobrut=hdfinfo([source '\' Dir_source(f).name]);
        dsetbrut=infobrut.SDS;
        brut=double(hdfread(dsetbrut));
        xx=find(brut(:)<-1);   %il faut conserver les "-1" pcq ils représentent les NaN
        brut(xx)=brut(xx)+256;
        pour = pourcentage_mod(brut);
       
        if pour>pour_min  
       
            %-------------------------%image des edge %----------------------------%
            % 0 = pixel non front / -1 = pixel front

            infosied=hdfinfo([source '\' Dir_source(f+1).name]);
            dsetsied=infosied.SDS;
            edge = hdfread(dsetsied);
            %------------------------------------------------------------------

            
            %CONVERSION NdG À DEGRE CELSIUS
            brut(brut==-1)=-4;
            %brut(brut==-1)=NaN;                              
            %brut(isnan(brut))=-4;
          
          
             brut(brut~=-4)=(brut(brut~=-4))*m+b;
              
          


            %--------------------------------------------------------------
            
            %SUPERPOSITION DES 2 IMAGES AVEC LE NETTOYAGE POUR LES COTES ET
            %LES NUAGES ET LA GLACE (FAUX FRONTS)
            dim_edge =size(edge );
            image_overlay=immerge_ATLANTIC(brut,edge ); % image_overlay remplace edge_complete de la version originale
            

            %-------FIGURE---------------------------------------
            image_overlay(image_overlay==-4)=land(image_overlay==-4)+image_overlay(image_overlay==-4);
            
            h=figure(1);
        
            %mapshow(lon,lat,image_overlay,'DisplayType','surface')  %       DECOMMENTER POUR PROJECTION
            imagesc(lon(1,:)',lat(:,1),image_overlay) %           COMMENTER POUR PROJECTION
            set(gca,'YDir','normal') %          COMMENTER POUR PROJECTION
            %hold all              %           DECOMMENTER POUR PROJECTION
            patch(coast(1:end-1,1),coast(1:end-1,2),[0.75 0.75 0.75])      %DECOMMENTER POUR ajouter le trait de cote
            colormap(mymap)
            %axis([ -80 -35 30 65]) %           DECOMMENTER POUR PROJECTION
            caxis([-2 30])
            colorbar
            title(['Overlay ',Dir_source(f).name(1:8)]) % nbre de caract total du nom du fichier de sortie
            
     
        % le répertoire traitement_OUTPUT\yyyy doit avoir été préalablement créé
            save([dest 'Overlay_' Dir_source(f).name(1:end-4) '.mat'],'image_overlay')
            saveas(h,[dest 'Overlay_' Dir_source(f).name(1:8)],'png')
            saveas(h,[dest 'Overlay_' Dir_source(f).name(1:8)],'fig')
            
            delete(h)
            
        end
    end
end

%         
