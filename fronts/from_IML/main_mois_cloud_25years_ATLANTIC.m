%   Guillaume Desbiens le 2 décembre 2011
%
%   Calcul la fréquence mensuelle de la couverture nuageuse, moyennée sur 1986-2010 
%   Les matrices sources se trouvent sous:J:\matlab_atlantic\moyenne_INPUT\
%   Les résulats sont enregistrés sous
%   J:\matlab_atlantic\probabilite_OUTPUT\cloud
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

clc
clear all
close all
dim_im=[3499 3000]; % taille de l'image 

%Pour mettre le fond blanc dans l'image
mymap=colormap(jet(50));
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


for m=1   % A MODIFIER selon les mois à traiter
    month=m;
    if  month==1
        strmonth2='January';
    elseif month==2    
        strmonth2='February';
    elseif month==3    
        strmonth2='March';
    elseif month==4    
        strmonth2='April';
    elseif month==5    
        strmonth2='May';
    elseif month==6    
        strmonth2='June';
    elseif month==7  
        strmonth2='July';
    elseif month==8   
        strmonth2='August';
    elseif month==9    
        strmonth2='September';
    elseif month==10    
        strmonth2='October';
    elseif month==11    
        strmonth2='November';
    elseif month==12    
        strmonth2='December';
    end    
    if (month<10)        
        strmonth= ['0' num2str(month)];
    else
        strmonth= num2str(month);
    end
   %------------------------------------------------------------------------
    source = 'J:\matlab_atlantic\moyenne_INPUT\';
   %------------------------

   % le fichier de destination DOIT ÊTRE CRÉÉ AVANT l'éxécution 

   dest ='probabilite_OUTPUT\cloud\';                                   
   
   [monthly_cloud_prob_25years,cumul_prob,used_files25] = prob_mois_cloud_25years_ATLANTIC(strmonth);
    

   for j=1:dim_im(1)
        for k=1:dim_im(2)
            if (land(j,k)== -1)
                monthly_cloud_prob_25years(j,k) = -1;          
            end
        end
   end

    h=figure(1);
    imagesc(lon(1,:)',lat(:,1),monthly_cloud_prob_25years)
    set(gca,'YDir','normal')

    patch(coast(1:end-1,1),coast(1:end-1,2),[0.75 0.75 0.75]) % ajoute le trait de côte à l'image

    caxis([-1 100])   % Echelle de couleur 0 à 100% 
    colormap(mymap)
    colorbar


    title([strmonth2 ' cloud frequency averaged over 1986-2010 (%) '  ]) 


    save([dest 'Cloud_' strmonth2 '_1986-2010.mat'],'lat','lon','used_files25','monthly_cloud_prob_25years','cumul_prob')                                                                            
    saveas(h,[dest 'Cloud_' strmonth2 '_1986-2010 '],'png') 
    saveas(h,[dest 'Cloud_' strmonth2 '_1986-2010 '],'fig')


    delete(h)
    clear('monthly_prob_25years','cumul_prob') % COMMENTER pour que les variables s'affichent 
                                               % (valider au moins un fichier, la variable used_files25 devrait être == 25)
                                           
end

