% Ce script va chercher les images recherchées selon l'intervalle de temps
% désiré (K:\Fronts\Edge_Results\yyyy) et créer une carte de probabilités de présence de fronts('prob_output\Mensuel\').

clc
clear all
close all
dim_im=[600 1100];

%Pour mettre le fond blanc dans l'image
mymap=colormap(jet);

mymap(1,:)=[0.9 0.9 0.9];
% mymap(2:4,:)=[1 1 1;1 1 1;1 1 1];


%-----Caractéristiques de la région gulf---------------
%-----LATITUDE/LONGITUDE------------------------------%

info = hdfinfo('beaufort.latlon.hdf');
dsets = info.SDS;
lat=hdfread(dsets(1));
lon=hdfread(dsets(2));
lat=(cutapq(lat,100,385,334,555))/100;
lon=cutapq(lon,100,385,334,555)/100;

%-----TRAIT DE COTE------------------------------%

coast=load('20728.dat');
infol=hdfinfo('landmask_beaufort_calypsopass.hdf');
dsetl=infol.SDS;
land = hdfread(dsetl);
land(land==1)=-1;
land=double(land);

%--------------------------------------------------------------------------
% Spécifier l'intervalle temporel de recherche des images
%--------------------------------------------------------------------------
 
for y=1998:2008  % 1998:2008 % A MODIFIER
    year=num2str(y)
    for m=10:10    % 5:10 % A MODIFIER
        month=m;
        if (month<10)        
            strmonth= ['0' num2str(month)];
        else
            strmonth= num2str(month);
        end
        %------------------------------------------------------------------------
        source = ['K:\Fronts\Edge_Results\' year  '\'];

        dest1 = ['prob_output\Mensuel\']; 
        %--------------------------------------------------------------------------
        % Aquisition des images appropriés dans les dossiers
        % K:\Fronts\Edge_Results\yyyy et creation de la carte probabilité de
        % présence de fronts
        %--------------------------------------------------------------------------

        mkdir(dest1,[year '_' strmonth])
        dest = [dest1 year '_' strmonth '\'];
       
        [c,count_pix,count_edge]=GetPicT(source,dest1,strmonth,dim_im); 

%         count_edge(count_edge==0)=land(count_edge==0)+count_edge(count_edge==0);  % pour mettre trait de cote
        
        h=figure(1);
        imagesc(lon(1,:)',lat(:,1),count_edge)
        set(gca,'YDir','normal')
        
%         patch(coast(1:end-1,1),coast(1:end-1,2),[0.75 0.75 0.75])

%        caxis([0 100])   % Echelle de couleur 0 a 100 a voir 
        colormap(mymap)
        colorbar
        title(['Probabilité de présence de front (%) ' dest(13:end) ' (' num2str(c) ')'])

        save([dest 'Prob_' dest(13:19) '_' year '_' strmonth '.mat'],'count_edge','count_pix','c')
        saveas(h,[dest 'Prob_' dest(13:19) '_' year '_' strmonth],'png')
        saveas(h,[dest 'Prob_' dest(13:19) '_' year '_' strmonth],'fig')
        
%         delete(h)
    end
end
