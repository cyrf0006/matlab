% ----------------------------
% Guillaume Desbiens 30 novembre 2011
%
% Ce script cr�er la matrice de probabilit� de pr�sence de fronts pour
% chacun des mois, moyenn�e sur 5 ans (2006-2010)

function [monthly_prob_5years,cumul_prob,used_files5] = prob_mois_5years_ATLANTIC(strmonth)

dim_im=[3499 3000];

source = 'J:\matlab_atlantic\moyenne_INPUT\';

Dir_source = dir(source); %utilise les .mat cr��s par "probabilites_MOIS_ATLANTIC.m" comme fichier d'entr�e
dim_source = size(Dir_source); 

cumul_prob = zeros(dim_im(1),dim_im(2));       % cr�ation de la matrice qui cumule les probabilit�s de d�tection de fronts de tous les mois, toutes les valeurs sont = 0
monthly_prob_5years=zeros(dim_im(1),dim_im(2)); % cr�ation de la matrice de probabilit�s de d�tection de fronts pour les 5 ans

extension  = 'mat';
used_files5=0; %compteurs de fichiers utilis�s; pour validation doit �tre = 5 (1 mois * 5 ans) 
     

for i=3:dim_source(1) 
    for y=2006:2010;
        year=num2str(y);
                  
            if (strcmp(Dir_source(i).name(14:17),year)&& strcmp(Dir_source(i).name(19:20),strmonth)&& strcmp(Dir_source(i).name(end-2:end),extension)) 
%                 s'assure que: - les caract�res 14 � 17 du nom du fichier source correspond � l'une des ann�es voulues
%                               - les caract�res 19 � 20 du nom du fichier source(Prob_Mensuel_1986_04.mat) correspondent au mois trait�
%                               - l'extension est .mat
               data=load([source Dir_source(i).name]);
               input=data.probability; 
               used_files5=used_files5+1;
               for j=1:dim_im(1)
                    for k=1:dim_im(2)
                        cumul_prob(j,k)=cumul_prob(j,k)+input(j,k);
                        monthly_prob_5years(j,k) = (cumul_prob(j,k)./used_files5);         
                    end
               end
            end                    
    end    
end
 



                