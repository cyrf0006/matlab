function [mean_prob_year,cumul_prob,used_files12] = prob_year_ATLANTIC(year)

% Ce script cr�er la matrice de: -  probabilit� de pr�sence de fronts selon l'ann�e voulue 
 


dim_im=[3499 3000];

%source = 'J:\matlab_atlantic\moyenne_INPUT\'; %  format =>> Prob_Mensuel_1990_02.mat
source= 'J:\matlab_atlantic\moyenne_INPUT\';
Dir_source = dir(source); %utilise les .mat cr��s par "probabilites_MOIS_ATLANTIC.m" comme fichier d'entr�e
dim_source = size(Dir_source); 

cumul_prob = zeros(dim_im(1),dim_im(2)); % cr�ation de la matrice qui cumule les prob de tous les mois, toutes les valeurs sont = 0
mean_prob_year=zeros(dim_im(1),dim_im(2)); % cr�ation de la matrice prob_season, toutes les valeurs sont = 0


extension  = 'mat';
used_files12=0; %compteurs de fichiers utilis�s; pour validation 

% 
for i=3:dim_source(1) 
    for m=1:12
        month=m;
          if (month<10)        
                strmonth= ['0' num2str(month)];
          else
                strmonth= num2str(month);
          end
                                    %name= Prob_Mensuel_1999_02.mat
          if (strcmp(Dir_source(i).name(14:17),year)&&strcmp(Dir_source(i).name(19:20),strmonth)&& strcmp(Dir_source(i).name(end-2:end),extension)) % strcmp : string comparaison
%                s'assure que: - les caract�res 14 � 17 du nom du fichier source correspond � l'une des ann�es voulues  
%                              - les caract�res 19 � 20 du nom du fichier source correspond � l'un des mois voulus
%                              - que l'extension est .mat
               data=load([source Dir_source(i).name]);
               input=data.probability;   
               used_files12=used_files12+1;

               for j=1:dim_im(1)
                        for k=1:dim_im(2)
                            cumul_prob(j,k)=cumul_prob(j,k)+input(j,k);
                            mean_prob_year(j,k) = (cumul_prob(j,k)./used_files12);
                        end
               end
           end
     end
 end
%             

    %---------CALCUL DE LA FR�QUENCE NUAGES/GLACES en %

    % for j=1:dim_im(1)
    %     for k=1:dim_im(2)
    %         if (count_pix(j,k)==0)
    %             cloud_prob(j,k) = 0;
    %         else cloud_prob(j,k)=((c-count_pix(j,k))/c)*100;   
    %         end
%      end
% end

                