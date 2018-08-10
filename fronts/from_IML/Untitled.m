clear all
close all
clc;

load J:\matlab_atlantic\probabilite_OUTPUT\mensuelle\'Prob_April_1986-2010.mat';



% imwrite('Prob_Mensuel_2002_08.mat.mat','Prob_Mensuel_2002_08.mat.hdf',hdf);
% imwrite(Prob_Mensuel_2002_08.mat,probability,'Prob_Mensuel_2002_08.mat.hdf',hdf);

imwrite(monthly_prob_25years,'Prob_April_1986-2010.hdf');


% imwrite(Prob_April_1986-2010,monthly_prob_25years,lat,lon,'your_hdf_file.hdf','Compression','none');