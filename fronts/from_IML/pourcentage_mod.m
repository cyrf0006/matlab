function pour = pourcentage_mod(image)
% Calcul le pourcentage d’eau sur l’image
dim = size(image);
good = 0; bad = 0;
for i=1:dim(1)
    for j=1:dim(2)
        if image(i,j)==-1
            bad = bad +1;
        else
            good = good +1;
        end
    end
end
pour = (good/(bad+good))*100;