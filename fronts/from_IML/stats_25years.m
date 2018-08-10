%%%%%%%%%%%%
% Ce script calcule la moyenne et l'écart-type des probabilités de chacune 
% des années traitées annuelle
close all;
clear all;
clc;
source='J:\matlab_atlantic\probabilite_OUTPUT\annuelle';

data=load('J:\matlab_atlantic\probabilite_OUTPUT\annuelle\Prob_1986-2010.mat');
data2=load('J:\matlab_atlantic\traitement_OUTPUT\2009\Overlay_20090730...L3..global.msst.mean.1d_2_Cut.mat');

mask=data2.image_overlay;
data25=data.mean_prob_25years;
data_25years=data25;

  
terre=find(mask(:)==-5); 
data_25years(terre)=data_25years(terre)-5; %attribut la valeur -5 aux pixels qui représentent la terre sur l'image 'image_overlay'

        
dimension=[3499 3000];
sum=0;
pix_inclus=0;
pix_exclus=0;

%calcul de la somme des pixels valides (différents de Terre)
for j=1:dimension(1)
    for k=1:dimension(2)
        if (data_25years(j,k)<0)
            pix_exclus = pix_exclus+1;
        else
             sum=sum+data_25years(j,k);
            
        end
    end
end


%calcul du nbre de pixels inclus ds les calculs
for j=1:dimension(1)
    for k=1:dimension(2)
        if (data_25years(j,k)>=0)
            pix_inclus = pix_inclus+1;
        end
    end
end

pix_checked=pix_inclus+pix_exclus; %calcul des pix visités
nb_pix=3499*3000;% nbre de pix dans l'image

total_mean=sum/(pix_inclus);% calcul de la moyenne de l'image 86-2010


% Calcul de l'écart type pour l'image 1986-2010
sum_sigma25_years=0;
sigma_pix=0;

for j=1:dimension(1)
    for k=1:dimension(2)
        if (data_25years(j,k)>=0)
            sigma_pix = (data_25years(j,k)-total_mean)^2; %calcul le sigma de chaque pixel
        end
        if (data_25years(j,k)>=0)
            sum_sigma25_years=sum_sigma25_years+sigma_pix; % calcul la somme des sigmas
        end
    end
end

sigma25_years=((1/pix_inclus)*sum_sigma25_years)^0.5; %clacul le sigma pour l'image entière (moyenne totale 1986-2010)








