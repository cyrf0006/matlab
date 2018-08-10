% ----------------------------
% Guillaume Desbiens 30 novembre 2011
%
% Ce script créer la matrice de probabilité de présence de fronts pour
% l'hiver, moyennée sur 5 ans (2006-2010)

function [mean_prob_season,cumul_prob,used_files15] = ECprob_winter_5years_ATLANTIC()

dim_im=[3499 3000];

source = 'J:\matlab_atlantic\moyenne_INPUT\'; 
Dir_source = dir(source); %utilise les .mat créés par "probabilites_MOIS_ATLANTIC.m" comme fichier d'entrée
dim_source = size(Dir_source); 

cumul_prob = zeros(dim_im(1),dim_im(2));       % création de la matrice qui cumule les probabilités de détection de fronts de tous les mois, toutes les valeurs sont = 0
mean_prob_season=zeros(dim_im(1),dim_im(2)); % création de la matrice de probabilités de détection de fronts pour la saison

extension  = 'mat';
used_files15=0; %compteurs de fichiers utilisés; pour validation doit être = à 15 (3 mois * 5 ans) 
      

for i=3:dim_source(1) 
    for y=2006:2010;
        year=num2str(y);
        for m=1:12
            month=m;
            if m==1||m==2||m==12
                if (month<10)        
                strmonth= ['0' num2str(month)];
                else
                    strmonth= num2str(month);
                end
                                        %name= Prob_Mensuel_1999_02.mat
                if (strcmp(Dir_source(i).name(14:17),year)&&strcmp(Dir_source(i).name(19:20),strmonth)&& strcmp(Dir_source(i).name(end-2:end),extension)) % strcmp : string comparaison
    %                s'assure que: - les caractères 14 à 17 du nom du fichier source correspond à l'un des mois voulus  
    %                              - les caractères 19 à 20 du nom du fichier source correspond à l'une des années voulues
    %                              - que l'extension est .mat

                   data=load([source Dir_source(i).name]);
                   input=data.probability; 
                   used_files15=used_files15+1;
                   for j=1:dim_im(1)
                        for k=1:dim_im(2)
                            cumul_prob(j,k)=cumul_prob(j,k)+input(j,k);
                            mean_prob_season(j,k) = (cumul_prob(j,k)./used_files15);       
                        end
                   end
                end  
            end    
        end    
    end
end 


                