% ----------------------------
% Guillaume Desbiens 30 novembre 2011
%
% Ce script cr�er la matrice de probabilit� de pr�sence de fronts pour
% l'hiver, moyenn�e sur 5 ans (2006-2010)

function [mean_prob_season,cumul_prob,used_files15] = ECprob_winter_5years_ATLANTIC()

dim_im=[3499 3000];

source = 'J:\matlab_atlantic\moyenne_INPUT\'; 
Dir_source = dir(source); %utilise les .mat cr��s par "probabilites_MOIS_ATLANTIC.m" comme fichier d'entr�e
dim_source = size(Dir_source); 

cumul_prob = zeros(dim_im(1),dim_im(2));       % cr�ation de la matrice qui cumule les probabilit�s de d�tection de fronts de tous les mois, toutes les valeurs sont = 0
mean_prob_season=zeros(dim_im(1),dim_im(2)); % cr�ation de la matrice de probabilit�s de d�tection de fronts pour la saison

extension  = 'mat';
used_files15=0; %compteurs de fichiers utilis�s; pour validation doit �tre = � 15 (3 mois * 5 ans) 
      

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
    %                s'assure que: - les caract�res 14 � 17 du nom du fichier source correspond � l'un des mois voulus  
    %                              - les caract�res 19 � 20 du nom du fichier source correspond � l'une des ann�es voulues
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


                