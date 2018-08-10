%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ce script calcule la moyenne et l'�cart-type des probabilit�s mensuelle
% pour chacune des ann�es
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

    dest ='J:\matlab_atlantic\probabilite_OUTPUT\stats\annuelle'; 

  