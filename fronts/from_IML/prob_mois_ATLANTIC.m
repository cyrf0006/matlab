function [image_overlay,c,count_pix,count_edge,probability,cloud_prob]=prob_mois_ATLANTIC(source,dest,strmonth,dim_im) 

% Ce script créer la matrice de: -  probabilité de présence de fronts 
%                                -  fréquence de nuage en % ex: le pix (i,j)était recouvert par des nuages 65% du temps, durant le mois xx
%
% selon l'intervalle de temps désiré.

%source = ['J:\matlab_atlantic\traitement_OUTPUT\' year  '\']; 
Dir_source = dir(source); %utilise les .mat créés par "traitement_edge" comme fichier d'entrée
                        %  source=['J:\matlab_atlantic\traitement_OUTPUT\'year '\']
dim_source = size(Dir_source); 
count_edge = zeros(dim_im(1),dim_im(2)); % création de la matrice count_edge, toutes les valeurs sont = 0
count_pix  = zeros(dim_im(1),dim_im(2));
probability= zeros(dim_im(1),dim_im(2));
cloud_prob = zeros(dim_im(1),dim_im(2));
extension  = 'mat';
c=0;

for i=3:dim_source(1)

        if (strcmp(Dir_source(i).name(13:14),strmonth)&& strcmp(Dir_source(i).name(end-2:end),extension)) % strcmp : string comparaison
           % s'assure que les 13 et 14 ième caractères du nom de la source correspond au mois && que l'extension est .mat
           data=load([source Dir_source(i).name]);
           image_overlay=data.image_overlay;   
            prob_input=data.image_overlay;
                I=(prob_input==-2.5);
                count_edge =count_edge+I;
            P=((prob_input~=-4)&(prob_input~=-5)); % ne compte que les pixels qui ne sont ni NaN(nuage/glace=> (-4)) ni de la terre=>(-5)
                count_pix=count_pix+P; 
            c=c+1;
        end
        P=0;
        I=0;
end

%----------CALCUL DE LA PROBABILITÉ-------------
for j=1:dim_im(1)
    for k=1:dim_im(2)
        if (count_pix(j,k)==0)
            count_edge(j,k) = 0;
        else
            if (count_edge(j,k)==0)
                count_edge(j,k) = 0;
            else
                probability(j,k) = (count_edge(j,k)./count_pix(j,k)).*100;

%                 if (probability(j,k)<=2)
%                     probability(j,k)=0;
%                 end 
            end
        end
    end
end
%----------------------------------------------

%---------CALCUL DE LA FRÉQUENCE NUAGES/GLACES en %

for j=1:dim_im(1)
    for k=1:dim_im(2)
        if (count_pix(j,k)==0)
            cloud_prob(j,k) = 0;
        else cloud_prob(j,k)=((c-count_pix(j,k))/c)*100;   
        end
     end
end

                