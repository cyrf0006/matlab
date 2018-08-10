%%%%%%%%%%%%
% Ce script calcule la moyenne et l'écart-type des probabilités mensuelle
close all;
clear all;
clc;

data=load('Prob_Mensuel_1990_12.mat');
data2=load('Overlay_20091230...L3..global.msst.mean.1d_2_Cut.mat');

mask=data2.image_overlay;
data06=data.probability;
data_clean06=data06;

  
terre=find(mask(:)==-5); 
data_clean06(terre)=data_clean06(terre)-5; %attribut la valeur -5 aux pixels qui représentent la terre sur l'image 'image_overlay'

        
dimension=[3499 3000];
sum=0;
pix_inclus=0;
pix_exclus=0;

%calcul de la somme des pixels valides (>0)
for j=1:dimension(1)
    for k=1:dimension(2)
        if (data_clean06(j,k)<0)
            pix_exclus = pix_exclus+1;
        else
             sum=sum+data_clean06(j,k);
            
        end
    end
end


%calcul du nbre de pixels inclus ds les calculs
for j=1:dimension(1)
    for k=1:dimension(2)
        if (data_clean06(j,k)>=0)
            pix_inclus = pix_inclus+1;
        end
    end
end

pix_checked=pix_inclus+pix_exclus; %calcul des pix visités
nb_pix=3499*3000;% nbre de pix dans l'image

moyenne06=sum/(pix_inclus);% calcul de la moyenne du mois

sum_sigma06=0;
sigma_pix=0;

for j=1:dimension(1)
    for k=1:dimension(2)
        if (data_clean06(j,k)>=0)
            sigma_pix = (data_clean06(j,k)-moyenne06)^2; %calcul le sigma de chaque pixel
        end
        if (data_clean06(j,k)>=0)
            sum_sigma06=sum_sigma06+sigma_pix; % calcul la somme des sigmas
        end
    end
end

sigma06=((1/pix_inclus)*sum_sigma06)^0.5; %clacul le sigma pour l'image entière (mensuel)








