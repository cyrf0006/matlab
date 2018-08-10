%anomalySST.m
% Calcul des anomalies de probabilité de détection de fronts:Ce programme va
% calculer les anomalies normalisées de probabilité de détection de fronts 
% mensuelles par rapport a la climatologie (1998-2008)et va tracer les cartes
% d'anomalie

clear all
close all
clc

year='2004'; %-----A MODIFIE ------
month='09'; %-----A MODIFIE ------
SD=0.774;    %-----A MODIFIE ------
%(SD10=0.796; SD09=0.774; SD08=1.25; SD07=1.778; SD06=0.898; SD05=0.114)

%Ouverture du fichier moyenne mensuelle totale (sur le nombre d'annees)
data1=load(['Moyenne_output\Moyennes_totales\im_moy_tot_' month '.mat']);
%lecture de la temperature moyenne totale de toutes les années
moy_tot=data1.im_moy; 
%Ouverture du fichier moyenne mensuelle sur un mois
data2=load(['C:\Program Files\MATLAB\R2007b\work\Audrey\Moyenne_output\Moyennes_mensuelles\im_moy_' year '_' month '.mat']);
%lecture de la temperature
moy_month=data2.im_moy;

%Calcul de l'anomalie de temperature normalisée
%anomaly=moy_month-moy_tot;
anomalySD=(moy_month-moy_tot)/SD;
anomalySD(isnan(anomalySD))=0;

%on rajoute le fichier lat_lon pour la carte
coord=load('SSTglobal1d_corners.txt');
lon=coord(1):coord(3):coord(2);%va calculer la grille de coordonnees avec le delta
lat=coord(4):coord(6):coord(5);

[north,a]=shaperead('northatlantic\northatlantic.shp');

%sauvegarde de la matrice   
save(['AnomalySD_output\anomalySD_' month '_'  year '.mat'],'anomalySD','lat','lon')

%creation de la carte et visualisation
figure;
imagesc(lon,lat,anomalySD)
set(gca,'YDir','normal')
hold on
mapshow(north,'Facecolor',[0.9 0.9 0.9])
caxis([-4 4])
colorbar
title(['SST monthly anomaly (Celsius): ' month '-' year '/(1998-2008)'])
xlabel('Longitude(°)')
ylabel('Latitude(°)')

%creation de la carte et visualisation
fh = figure(2);
imagesc(lon,lat,anomalySD)
set(gca,'YDir','normal')
hold on
mapshow(north,'Facecolor',[0.9 0.9 0.9])
caxis([-4 4])
title(['SST monthly anomaly (Celsius): ' month '-' year '/(1998-2008)'])
xlabel('Longitude(°)')
ylabel('Latitude(°)')
load('mymap.txt')
colormap(mymap)
colorbar

%sauvegarde en epsc : meilleur format qui garde le fond blanc et les couleurs de l'image

saveas(fh,['AnomalySD_output\anomaly_' month '_' year],'epsc')
saveas(fh,['AnomalySD_output\anomaly_' month '_' year],'png')
saveas(fh,['AnomalySD_output\anomaly_' month '_' year],'fig')

% %on peut aussi sauvegarder en fig et apres on ouvre le fig et on va faire
% %dans Edit: copy figure et la mettre directement dans un ppt par exemple(c'est une alternative)
% min(min(anomaly))
% max(max(anomaly))
% 
% %-----------------------------------------------------------------
%        %----------------------------------------------------------
%                    %----------------------------------------------
%                    %----------------------------------------------
%                    %ou bien appliquer ce programme pour avoir la projection conique 
% 
% %nom du fichier contenant la matrice de latitude et longitude
% nom2=['C:\Program Files\MATLAB\R2007b\work\beaufort.latlon.hdf'];
% %mise a zero
% scale_latlong=[];latitude=[];longitude=[];
% %lecture de l'échelle, la latitude et longitude(avec application de l'échelle)
% scale_latlong=hdfread(nom2,'scale_factor');
% latitude=double(hdfread(nom2,'latitude'))*scale_latlong{1,1};
% longitude=double(hdfread(nom2,'longitude'))*scale_latlong{1,1};
% %ouverture de la figure
% carte=[];
% h=figure;
% %arrangement des axes
% objet1=[];
% objet1=axes('DataAspectRatioMode','manual','DataAspectRatio',[0.5 1 1]);
% mapshow(lon(1:end-1),lat(1:end-1),anomaly,'DisplayType','surface') %changer le nom en anomaly % "surface" il save en png et texturemap save en fig
% hold all
% grid on
% %affichage de la barre de couleur
% objet4=[];
% objet4=colorbar;
% 
% mymap=colormap('jet');
% mymap(1,:)=[0.75 0.75 0.75];
% colormap(mymap)
% 
% xlabel('Longitude(°)')
% ylabel('Latitude(°)')
% 
% %titre
% objet5=[];
% objet5=title(['Anomaly SST 07/1998']);
% %on arrange les dimentions
% set(objet1,'DataAspectRatio',[1 0.5 1])
% %sauvegarde de la carte
% %print(h,'-djpeg','-r500',['anomalycubic07_1998.png'])
% 
