% ------------------------------
% Ce programme permet de visualiser les valeurs des pixels des images .hdf en CN
%--------------------------------------------------------------------------


% landmask_hudson_binary.hdf est le fichier à visualiser
% Binarize est le 'name' du dataset, qui est contenu dans la variable info_SDS


clc 
clear all
close all

info=hdfinfo('19860707...L3..global.msst.mean.1d_3_Cut.hdf')
info_SDS=info.SDS
dset=info.SDS
input_visual=hdfread('19860707...L3..global.msst.mean.1d_3_Cut.hdf','<CutImg>');

% input_visual = int8(-100);
% b = cast(input_visual,'uint8');
% class(b)






% right-click/imagesc(data) sur la variable 'data' permet de voir l'image 
% et ensuite avec l'outil 'data cursor' d'obtenir la valeur en CN (index)
% en double-clickant sur data on peut consulter la matrice