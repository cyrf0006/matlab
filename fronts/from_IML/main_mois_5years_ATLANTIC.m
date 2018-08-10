%   Guillaume Desbiens le 2 d�cembre 2011
%
%   Calcul des moyennes des probabilit�s de presence de fronts
%   Ce programme va chercher les matrices de probabilit�s mensuelles de d�tection de 
%   fronts et va calculer les climatologies selon le mois choisi, et ce, sur la p�riode 2006-2010. 
%   Les matrices sources se trouvent sous:J:\matlab_atlantic\moyenne_INPUT\
%   Les r�sulats sont enregistr�s sous
%   J:\matlab_atlantic\probabilite_OUTPUT\mensuelle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

clc
clear all
close all
dim_im=[3499 3000]; % taille de l'image 

%Pour mettre le fond blanc dans l'image
mymap=colormap(jet(50));
mymap(1,:)=[0.9 0.9 0.9];
% mymap(2:4,:)=[1 1 1;1 1 1;1 1 1];


%-----Caract�ristiques de la r�gion gulf---------------
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


for m=1:12   % A MODIFIER selon les mois � traiter
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

   % le fichier de destination DOIT �TRE CR�� AVANT l'�x�cution 

   dest ='probabilite_OUTPUT\mensuelle\';                                   
   
   [monthly_prob_5years,cumul_prob,used_files5] = prob_mois_5years_ATLANTIC(strmonth);
    

    h=figure(1);
    imagesc(lon(1,:)',lat(:,1),monthly_prob_5years)
    set(gca,'YDir','normal')

    patch(coast(1:end-1,1),coast(1:end-1,2),[0.75 0.75 0.75]) % ajoute le trait de c�te � l'image

    caxis([0 20])   % Echelle de couleur 0 a 20% de probabilit� (les valeurs >20 sont toutes de la m�me couleur)
    colormap(mymap)
    colorbar


    title([strmonth2 ' front detection probability averaged over 2006-2010 (%) '  ]) 


    save([dest 'Prob_' strmonth2 '_2006-2010.mat'],'lat','lon','used_files5','monthly_prob_5years','cumul_prob')                                                                            
    saveas(h,[dest 'Prob_' strmonth2 '_2006-2010 '],'png') 
    saveas(h,[dest 'Prob_' strmonth2 '_2006-2010 '],'fig')
    
   
% POUR UNE RAISON OBSCURE SI LA PARTIE CR�ATION DES FICHIERS HDF EST EX�CUT�E, LE MESSAGE
% D'ERREUR SUIVANT APPARA�T (pour la ligne 82):
    
%??? Undefined function or method 'prob_mois_25years_ATLANTIC' for input arguments of type 'char'.
   
%    prob_month5years_integer=uint8(monthly_prob_5years);% transforme la matrice monthly_prob_5years en integer (4.4 devient 4)
%                                                          % n�cessaire pour le fichier HDF
% 
%     %-----  Cr�ation des fichiers HDF
%     
%     cd('J:\matlab_atlantic\probabilite_OUTPUT\HDF\monthly\');  
% 
%     if  month==1
%             imwrite(prob_month5years_integer,'Prob_January_2006-2010.hdf');
%         elseif month==2    
%             imwrite(prob_month5years_integer,'Prob_February_2006-2010.hdf');
%         elseif month==3    
%              imwrite(prob_month5years_integer,'Prob_March_2006-2010.hdf');
%         elseif month==4    
%              imwrite(prob_month5years_integer,'Prob_April_2006-2010.hdf');
%         elseif month==5    
%              imwrite(prob_month5years_integer,'Prob_May_2006-2010.hdf');
%         elseif month==6    
%              imwrite(prob_month5years_integer,'Prob_June_2006-2010.hdf');
%         elseif month==7  
%              imwrite(prob_month5years_integer,'Prob_July_2006-2010.hdf');
%         elseif month==8   
%              imwrite(prob_month5years_integer,'Prob_August_2006-2010.hdf');
%         elseif month==9    
%              imwrite(prob_month5years_integer,'Prob_September_2006-2010.hdf');
%         elseif month==10    
%              imwrite(prob_month5years_integer,'Prob_October_2006-2010.hdf');
%         elseif month==11    
%              imwrite(prob_month5years_integer,'Prob_November_2006-2010.hdf');
%         elseif month==12    
%              imwrite(prob_month5years_integer,'Prob_December_2006-2010.hdf');
%     end    

    delete(h)
    clear('monthly_prob_5years','cumul_prob') % COMMENTER pour que les variables s'affichent 
                                               % (valider au moins un fichier, la variable used_files5 devrait �tre == 25
end

