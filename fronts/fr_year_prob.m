function [image_overlay,c,count_pix,count_edge,probability,cloud_prob]=fr_year_prob(source,dest,strmonth,dim_im) 

% Ce script créer la matrice de: -  probabilité de présence de fronts 
%                                -  fréquence de nuage en % ex: le pix (i,j)était recouvert par des nuages 65% du temps, durant le mois xx
%
% selon l'intervalle de temps désiré.

%source = ['J:\matlab_atlantic\traitement_OUTPUT\' year  '\']; 

% Load 1st image to get size
load(sstFiles(1,:)) % load 1st image
    
% Matrices initalization    
edgeCount = zeros(=size(image_overlay)); % fronts detection
pixelCount  = zeros(=size(image_overlay)); % cloud/ice free pixel
probability= zeros(=size(image_overlay)); % probability
cloud_prob = zeros(=size(image_overlay)); % cloub prob.
imageCount = size(sstFiles, 1);

for i = 1:size(sstFiles, 1)    
    data = load(sstFiles(i,:));
    image_overlay=data.image_overlay;   
    prob_input=data.image_overlay;

    I=(prob_input==-2.5);
    edgeCount = edgeCount + I;
    
    P=((prob_input~=-4)&(prob_input~=-5)); % Ignore ICE (-4) and CONTINENT (-5)
    pixelCount = pixelCount + P; 
end

keyboard

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

                