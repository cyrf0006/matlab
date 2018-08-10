function [cloud_season,used_files2] = cloud_season_ATLANTIC()


% Ce script créer la matrice de: -  probabilité de présence de fronts selon la saIson voulue 
 
season=1; 

dim_im=[3499 3000];

%source = 'J:\matlab_atlantic\moyenne_INPUT\'; %  format =>> Prob_Mensuel_1990_02.mat
source= 'J:\matlab_atlantic\mmm\';
Dir_source = dir(source); %utilise les .mat créés par "probabilites_MOIS_ATLANTIC.m" comme fichier d'entrée
dim_source = size(Dir_source); 

% total_prob = zeros(dim_im(1),dim_im(2));       % création de la matrice qui cumule les probabilités de détection de fronts de tous les mois, toutes les valeurs sont = 0
% probability_season=zeros(dim_im(1),dim_im(2)); % création de la matrice de probabilités de détection de fronts pour la saison
total_cloud=zeros(dim_im(1),dim_im(2));        % création de la matrice qui cumule les probabilités de nuages de tous les mois,
cloud_season=zeros(dim_im(1),dim_im(2));       % création de la matrice de probabilités de nuage pour la saison

extension  = 'mat';
used_files2=0; %compteurs de fichiers utilisés; pour validation 

% if (season==1)        
%     month=1:3;
%     elseif (season==2)
%     month=4:6;
%     elseif (season==3)
%     month=7:9;
%     elseif strseason==4
%     month=10:12;
%         

% for y=1999  % A MODIFIER
% for season=1
    y=1999;
     year=num2str(y);
        for i=3:dim_source(1) 
            if (season==1)        
                for m=1:3
                    month=m;
                    if (month<10)        
                        strmonth= ['0' num2str(month)];
                    else
                        strmonth= num2str(month);
                    end
                                            %name= Prob_Mensuel_1999_02.mat
                    if (strcmp(Dir_source(i).name(14:17),year)&&strcmp(Dir_source(i).name(19:20),strmonth)&& strcmp(Dir_source(i).name(end-2:end),extension)) % strcmp : string comparaison
                       % s'assure que les 19 et 20 ième caractères du nom de la source correspond au mois && que l'extension est .mat
                       data=load([source Dir_source(i).name]);
                       cloud=data.cloud_prob;
                       used_files2=used_files2+1;
                           for j=1:dim_im(1)
                                for k=1:dim_im(2)
%                                     total_prob(j,k)=total_prob(j,k)+input(j,k);
%                                     probability_season(j,k) = (total_prob(j,k)./used_files);
                                    total_cloud(j,k)=total_cloud(j,k)+cloud(j,k);
                                    cloud_season(j,k)=(total_cloud(j,k)./used_files2);
                                end
                           end
                     end
                 end
            end
        end             
%             if (season==2)        
%                 for m=4:6
%                 %m=1:3;
%                 month=m;
%                     if (month<10)        
%                         strmonth= ['0' num2str(month)];
%                     else
%                         strmonth= num2str(month);
%                     end
%                                             %name= Prob_Mensuel_1999_02.mat
%                     if (strcmp(Dir_source(i).name(14:17),year)&&strcmp(Dir_source(i).name(19:20),strmonth)&& strcmp(Dir_source(i).name(end-2:end),extension)) % strcmp : string comparaison
%                        % s'assure que les 19 et 20 ième caractères du nom de la source correspond au mois && que l'extension est .mat
%                        data=load([source Dir_source(i).name]);
%                        input=data.probability;   
%                        used_files=used_files+1;
% 
%                        for j=1:dim_im(1)
%                             for k=1:dim_im(2)
%                                 total_prob(j,k)=total_prob(j,k)+input(j,k);
%                                 
%                             end
%                        end
%                     end
%                 end
%             end
% end
% probability_season(j,k) = (total_prob(j,k)./used_files);

                