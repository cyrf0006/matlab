% function[STD_sigma_year,sigma_year,sigma_25years,cumul_sigma,used_files25]=anomaly_year_ATLANTIC(source,dim_im,extension)
function [STD_sigma_year,sigma_year,used_files25]=anomaly_year_ATLANTIC()

source='J:\matlab_atlantic\probabilite_OUTPUT\annuelle\';
Dir_source = dir(source); 
dim_source = size(Dir_source); 
dim_im=[3499 3000];
%Ouverture du fichier moyenne annuelle totale 
data=load('J:\matlab_atlantic\probabilite_OUTPUT\annuelle\Prob_1986-2010.mat');

%lecture de la probabilit� moyenne totale de toutes les ann�es
mean_25y=data.mean_prob_25years;


% cr�ation de la matrice qui contient l'�cart-type de chaque pixel pour l'ann�e trait�e 
sigma_year=zeros(dim_im(1),dim_im(2));

% cr�ation de la matrice qui contient le cumulatif des �cart-type de chaque pixel pour 1986-2010  
cumul_sigma=zeros(dim_im(1),dim_im(2));

% cr�ation de la matrice qui contient l'�cart-type 1986-2010 de chaque pixel   
sigma_25years=zeros(dim_im(1),dim_im(2));

% cr�ation de la matrice qui contient l'�cart-type STANDARDIS� de chaque pixel pour l'ann�e trait�e 
STD_sigma_year=zeros(dim_im(1),dim_im(2));

used_files25=0;
extension='.mat';



for i=3:dim_source(1) 
    for y=1986:2010
        year=num2str(y);
                                    %name= Prob_1999.mat
            if (strcmp(Dir_source(i).name(6:9),year))&& (strcmp(Dir_source(i).name(end-2:end),extension)) % strcmp : string comparaison
%                s'assure que: - les caract�res 6 � 9 du nom du fichier source correspond � l'une des ann�es voulues  
%                              - que l'extension est .mat

%     PROBLEME A CE NIVEAU : LE USED_FILES25 RESTE � Z�RO !!

               data=load([source Dir_source(i).name]);
               mean_year=data.mean_prob_year; %charge la variable de la prob moyenne annuelle
               used_files25=used_files25+1; %compteur de fichiers utilis�s, devrait �tre ==0
               
                for j=1:dim_im(1)
                    for k=1:dim_im(2)
                            sigma_year(j,k)=mean_year(j,k)-mean_25y(j,k); % calcul de l'�cart-type de chaque pixel pour l'ann�e en cours de traitement 
                            cumul_sigma(j,k)=cumul_sigma(j,k)+sigma_year(j,k); % cumule les �-t de toutes les ann�es 
                    end
                end

                % calcul de l'�-t 1986-2010 pour chaque pixel         
                for j=1:dim_im(1)
                    for k=1:dim_im(2)
                        sigma_25years(j,k) =((1/used_files25(j,k))*cumul_sigma(j,k))^0.5; 
                    end
                end 

                % calcul de l'�-t STANDARDIS� pour chaque pixel de l'ann�e en cours de traitement  
                for j=1:dim_im(1)
                    for k=1:dim_im(2)
                        STD_sigma_year(j,k)=(mean_year(j,k)-mean_25y(j,k))/sigma_25years(j,k);
                    end
                end 
            end
            
            
     end           
end



