function[sigma_25years,cumul_sigma,used_files25]=stats_25years_ATLANTIC(source,dim_im,extension)


% source='E:\matlab_atlantic\probabilite_OUTPUT\annuelle';

%Ouverture du fichier moyenne annuelle totale 
data=load('E:\matlab_atlantic\probabilite_OUTPUT\annuelle\Prob_1986-2010.mat');
%lecture de la probabilit� moyenne totale de toutes les ann�es
mean_25y=data.mean_prob_25years

% data2=load('E:\matlab_atlantic\traitement_OUTPUT\2009\Overlay_20090730...L3..global.msst.mean.1d_2_Cut.mat');


% mask=data2.image_overlay;
% data25=data.mean_prob_25years;
% data_25years=data25;
% 
%   
% terre=find(mask(:)==-5); 
% data_25years(terre)=data_25years(terre)-5; %attribut la valeur -5 aux pixels qui repr�sentent la terre sur l'image 'image_overlay'

% cr�ation de la matrice qui contient l'�cart-type de chaque pixel 
sigma_year=zeros(dim_im(1),dim_im(2));
% cr�ation de la matrice qui contient le cumulatif des �cart-type de chaque pixel pour 1986-2010  
cumul_sigma=zeros(dim_im(1),dim_im(2));
% cr�ation de la matrice qui contient l'�cart-type 1986-2010 de chaque pixel   
sigma_25years=zeros(dim_im(1),dim_im(2));


used_files25=0
        
% dim_im=[3499 3000];
% sum=0;
pix_inclus=0;
pix_exclus=0;
% extension='.mat';

for i=3:dim_source(1) 
    for y=2009;
        year=num2str(y);
                                    %name= Prob_Mensuel_1999_02.mat
            if (strcmp(Dir_source(i).name(14:17),year))&& strcmp(Dir_source(i).name(end-2:end),extension)) % strcmp : string comparaison
%                s'assure que: - les caract�res 14 � 17 du nom du fichier source correspond � l'une des ann�es voulues  
%                              - que l'extension est .mat

               data=load([source Dir_source(i).name]);
               mean_year=data.mean_prob_year;
               used_files25=used_files25+1;
               
            for j=1:dim_im(1)
                for k=1:dim_im(2)
                        sigma_year=mean_year - mean_25y; % calcul de l'�cart-type de chaque pixel pour l'ann�e en cours de traitement 
                        cumul_sigma=cumul_sigma+sigma_year; % cumule les �-t de toutes les ann�es 
                end
            end
            
for j=1:dim_im(1)
    for k=1:dim_im(2)
            sigma_25years =((1/used_files25)*cumul_sigma)^0.5; 
    end
            end            
           

pix_checked=pix_inclus+pix_exclus; %calcul des pix visit�s
nb_pix=3499*3000;% nbre de pix dans l'image

total_mean=sum/(pix_inclus);% calcul de la moyenne de l'image 86-2010


% Calcul de l'�cart type pour l'image 1986-2010
sum_sigma25_years=0;
sigma_pix=0;

% for j=1:dimension(1)
%     for k=1:dimension(2)
%         if (data_25years(j,k)>=0)
%             sigma_pix = (data_25years(j,k)-total_mean)^2; %calcul le sigma de chaque pixel
%         end
%         if (data_25years(j,k)>=0)
%             sum_sigma25_years=sum_sigma25_years+sigma_pix; % calcul la somme des sigmas
%         end
%     end
% end

sigma25_years=((1/pix_inclus)*sum_sigma25_years)^0.5; %clacul le sigma pour l'image enti�re (moyenne totale 1986-2010)




