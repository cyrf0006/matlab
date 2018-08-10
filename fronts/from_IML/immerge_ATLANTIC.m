function [mpo_brut] = immerge_ATLANTIC(mpo_brut,mpo_edge_clean)
% Superpose les edge  sur l'image SST en mettant les edge  � -2.5 et les
% donn�es NaN (Terre, Nuage et Glace) � -4
dim_edge =size(mpo_edge_clean);

for i=2:(dim_edge (1)-1)
    for j=2:(dim_edge (2)-1)
        if(mpo_edge_clean(i,j)==-1)
                if (mpo_brut(i,j) == -4)
                elseif (mpo_brut(i+1,j) == -4)
                elseif (mpo_brut(i,j+1) == -4)
                elseif (mpo_brut(i-1,j) == -4)
                elseif (mpo_brut(i,j-1) == -4)
                elseif (mpo_brut(i+1,j+1) == -4)
                elseif (mpo_brut(i-1,j+1) == -4)
                elseif (mpo_brut(i+1,j-1) == -4)
                elseif (mpo_brut(i-1,j-1) == -4)
                else
                     mpo_brut(i,j)=-2.5;
                end
        end
    end
end