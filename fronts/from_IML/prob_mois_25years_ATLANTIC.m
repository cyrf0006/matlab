% ----------------------------
% Guillaume Desbiens 30 novembre 2011
%
% Ce script créer la matrice de probabilité de présence de fronts pour
% chacun des mois, moyennée sur 25 ans (1986-2010)

function [monthly_prob_25years,cumul_prob,used_files25] = prob_mois_25years_ATLANTIC(strmonth)

dim_im=[3499 3000];

source = 'J:\matlab_atlantic\moyenne_INPUT\';

Dir_source = dir(source); %utilise les .mat créés par "probabilites_MOIS_ATLANTIC.m" comme fichier d'entrée
dim_source = size(Dir_source); 

cumul_prob = zeros(dim_im(1),dim_im(2));       % création de la matrice qui cumule les probabilités de détection de fronts de tous les mois, toutes les valeurs sont = 0
monthly_prob_25years=zeros(dim_im(1),dim_im(2)); % création de la matrice de probabilités de détection de fronts pour les 25 ans

extension  = 'mat';
used_files25=0; %compteurs de fichiers utilisés; pour validation doit être = 25 (1 mois * 25 ans) 
     

for i=3:dim_source(1) 
    for y=1986:2010;
        year=num2str(y);
                  
            if (strcmp(Dir_source(i).name(14:17),year)&& strcmp(Dir_source(i).name(19:20),strmonth)&& strcmp(Dir_source(i).name(end-2:end),extension)) 
%                 s'assure que: - les caractères 14 à 17 du nom du fichier source correspond à l'une des années voulues
%                               - les caractères 19 à 20 du nom du fichier source(Prob_Mensuel_1986_04.mat) correspondent au mois traité
%                               - l'extension est .mat
               data=load([source Dir_source(i).name]);
               input=data.probability; 
               used_files25=used_files25+1;
               for j=1:dim_im(1)
                    for k=1:dim_im(2)
                        cumul_prob(j,k)=cumul_prob(j,k)+input(j,k);
                        monthly_prob_25years(j,k) = (cumul_prob(j,k)./used_files25);         
                    end
               end
            end                    
    end    
end
 



                