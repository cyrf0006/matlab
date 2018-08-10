clear all;
close all;
clc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  MOIS SUR 25 ANS          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load J:\matlab_hudson\probabilite_OUTPUT\mensuelle\'Prob_January_1986-2010.mat';
% 
% prob_8int=uint8(monthly_prob_25years);
% 
% cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\monthly\');   
% imwrite(prob_8int,'Prob_January_1986-2010.hdf');
% 
% %--------------------------
% 
%                      
% load J:\matlab_hudson\probabilite_OUTPUT\mensuelle\'Prob_February_1986-2010.mat';
% 
% prob_8int=uint8(monthly_prob_25years);
% 
% cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\monthly\');   
% imwrite(prob_8int,'Prob_February_1986-2010.hdf');
% 
% %--------------------------
% %--------------------------
% 
%                      
% load J:\matlab_hudson\probabilite_OUTPUT\mensuelle\'Prob_March_1986-2010.mat';
% 
% prob_8int=uint8(monthly_prob_25years);
% 
% cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\monthly\');   
% imwrite(prob_8int,'Prob_March_1986-2010.hdf');
% 
% %--------------------------%--------------------------
% 
%                      
% load J:\matlab_hudson\probabilite_OUTPUT\mensuelle\'Prob_April_1986-2010.mat';
% 
% prob_8int=uint8(monthly_prob_25years);
% 
% cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\monthly\');   
% imwrite(prob_8int,'Prob_April_1986-2010.hdf');

%--------------------------%--------------------------

                     
load J:\matlab_hudson\probabilite_OUTPUT\mensuelle\'Prob_May_1986-2010.mat';

prob_8int=uint8(monthly_prob_25years);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\monthly\');   
imwrite(prob_8int,'Prob_May_1986-2010.hdf');

%--------------------------
load J:\matlab_hudson\probabilite_OUTPUT\mensuelle\'Prob_June_1986-2010.mat';

prob_8int=uint8(monthly_prob_25years);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\monthly\');   
imwrite(prob_8int,'Prob_June_1986-2010.hdf');

%--------------------------
load J:\matlab_hudson\probabilite_OUTPUT\mensuelle\'Prob_July_1986-2010.mat';

prob_8int=uint8(monthly_prob_25years);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\monthly\');   
imwrite(prob_8int,'Prob_July_1986-2010.hdf');

%--------------------------
load J:\matlab_hudson\probabilite_OUTPUT\mensuelle\'Prob_August_1986-2010.mat';

prob_8int=uint8(monthly_prob_25years);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\monthly\');   
imwrite(prob_8int,'Prob_August_1986-2010.hdf');

%--------------------------
load J:\matlab_hudson\probabilite_OUTPUT\mensuelle\'Prob_September_1986-2010.mat';

prob_8int=uint8(monthly_prob_25years);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\monthly\');   
imwrite(prob_8int,'Prob_September_1986-2010.hdf');
%--------------------------
load J:\matlab_hudson\probabilite_OUTPUT\mensuelle\'Prob_October_1986-2010.mat';

prob_8int=uint8(monthly_prob_25years);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\monthly\');   
imwrite(prob_8int,'Prob_October_1986-2010.hdf');

%--------------------------
load J:\matlab_hudson\probabilite_OUTPUT\mensuelle\'Prob_November_1986-2010.mat';

prob_8int=uint8(monthly_prob_25years);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\monthly\');   
imwrite(prob_8int,'Prob_November_1986-2010.hdf');

%--------------------------
% load J:\matlab_hudson\probabilite_OUTPUT\mensuelle\'Prob_December_1986-2010.mat';
% 
% prob_8int=uint8(monthly_prob_25years);
% 
% cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\monthly\');   
% imwrite(prob_8int,'Prob_December_1986-2010.hdf');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  SAISONS         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load J:\matlab_hudson\probabilite_OUTPUT\saisonniere\'Prob_winter_1986-2010.mat';

prob_8int=uint8(season_prob_25years);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\monthly\');   
imwrite(prob_8int,'Prob_November_1986-2010.hdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  ANNEES          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load J:\matlab_hudson\probabilite_OUTPUT\annuelle\'Prob_1986.mat';
prob_8int=uint8(mean_prob_year);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\yearly\');   
imwrite(prob_8int,'Prob_1986.hdf');

%---------------------------------------
load J:\matlab_hudson\probabilite_OUTPUT\annuelle\'Prob_1987.mat';
prob_8int=uint8(mean_prob_year);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\yearly\');   
imwrite(prob_8int,'Prob_1987.hdf');

%---------------------------------------
load J:\matlab_hudson\probabilite_OUTPUT\annuelle\'Prob_1988.mat';
prob_8int=uint8(mean_prob_year);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\yearly\');   
imwrite(prob_8int,'Prob_1988.hdf');

%---------------------------------------
load J:\matlab_hudson\probabilite_OUTPUT\annuelle\'Prob_1989.mat';
prob_8int=uint8(mean_prob_year);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\yearly\');   
imwrite(prob_8int,'Prob_1989.hdf');

%---------------------------------------
load J:\matlab_hudson\probabilite_OUTPUT\annuelle\'Prob_1990.mat';
prob_8int=uint8(mean_prob_year);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\yearly\');   
imwrite(prob_8int,'Prob_1990.hdf');
%---------------------------------------
load J:\matlab_hudson\probabilite_OUTPUT\annuelle\'Prob_1991.mat';
prob_8int=uint8(mean_prob_year);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\yearly\');   
imwrite(prob_8int,'Prob_1991.hdf');
%---------------------------------------
load J:\matlab_hudson\probabilite_OUTPUT\annuelle\'Prob_1992.mat';

prob_8int=uint8(mean_prob_year);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\yearly\');   
imwrite(prob_8int,'Prob_1992.hdf');


%---------------------------------------
load J:\matlab_hudson\probabilite_OUTPUT\annuelle\'Prob_1993.mat';
prob_8int=uint8(mean_prob_year);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\yearly\');   
imwrite(prob_8int,'Prob_1993.hdf');


%---------------------------------------
load J:\matlab_hudson\probabilite_OUTPUT\annuelle\'Prob_1994.mat';

prob_8int=uint8(mean_prob_year);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\yearly\');   
imwrite(prob_8int,'Prob_1994.hdf');


%---------------------------------------
load J:\matlab_hudson\probabilite_OUTPUT\annuelle\'Prob_1995.mat';
prob_8int=uint8(mean_prob_year);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\yearly\');   
imwrite(prob_8int,'Prob_1995.hdf');


%---------------------------------------
load J:\matlab_hudson\probabilite_OUTPUT\annuelle\'Prob_1996.mat';

prob_8int=uint8(mean_prob_year);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\yearly\');   
imwrite(prob_8int,'Prob_1996.hdf');


%---------------------------------------
load J:\matlab_hudson\probabilite_OUTPUT\annuelle\'Prob_1997.mat';
prob_8int=uint8(mean_prob_year);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\yearly\');   
imwrite(prob_8int,'Prob_1997.hdf');


%---------------------------------------
load J:\matlab_hudson\probabilite_OUTPUT\annuelle\'Prob_1998.mat';
prob_8int=uint8(mean_prob_year);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\yearly\');   
imwrite(prob_8int,'Prob_1998.hdf');


%---------------------------------------
load J:\matlab_hudson\probabilite_OUTPUT\annuelle\'Prob_1999.mat';
prob_8int=uint8(mean_prob_year);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\yearly\');   
imwrite(prob_8int,'Prob_1999.hdf');


%---------------------------------------
load J:\matlab_hudson\probabilite_OUTPUT\annuelle\'Prob_2000.mat';
prob_8int=uint8(mean_prob_year);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\yearly\');   
imwrite(prob_8int,'Prob_2000.hdf');


%---------------------------------------
load J:\matlab_hudson\probabilite_OUTPUT\annuelle\'Prob_2001.mat';
prob_8int=uint8(mean_prob_year);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\yearly\');   
imwrite(prob_8int,'Prob_2001.hdf');


%---------------------------------------
load J:\matlab_hudson\probabilite_OUTPUT\annuelle\'Prob_2002.mat';
prob_8int=uint8(mean_prob_year);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\yearly\');   
imwrite(prob_8int,'Prob_2002.hdf');


%---------------------------------------
load J:\matlab_hudson\probabilite_OUTPUT\annuelle\'Prob_2003.mat';
prob_8int=uint8(mean_prob_year);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\yearly\');   
imwrite(prob_8int,'Prob_2003.hdf');


%---------------------------------------
load J:\matlab_hudson\probabilite_OUTPUT\annuelle\'Prob_2004.mat';
prob_8int=uint8(mean_prob_year);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\yearly\');   
imwrite(prob_8int,'Prob_2004.hdf');


%---------------------------------------
load J:\matlab_hudson\probabilite_OUTPUT\annuelle\'Prob_2005.mat';
prob_8int=uint8(mean_prob_year);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\yearly\');   
imwrite(prob_8int,'Prob_2005.hdf');
%---------------------------------------
load J:\matlab_hudson\probabilite_OUTPUT\annuelle\'Prob_2006.mat';
prob_8int=uint8(mean_prob_year);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\yearly\');   
imwrite(prob_8int,'Prob_2006.hdf');


%---------------------------------------
load J:\matlab_hudson\probabilite_OUTPUT\annuelle\'Prob_2007.mat';
prob_8int=uint8(mean_prob_year);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\yearly\');   
imwrite(prob_8int,'Prob_2007.hdf');
%---------------------------------------
load J:\matlab_hudson\probabilite_OUTPUT\annuelle\'Prob_2008.mat';
prob_8int=uint8(mean_prob_year);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\yearly\');   
imwrite(prob_8int,'Prob_2008.hdf');


%---------------------------------------
load J:\matlab_hudson\probabilite_OUTPUT\annuelle\'Prob_2009.mat';
prob_8int=uint8(mean_prob_year);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\yearly\');   
imwrite(prob_8int,'Prob_2009.hdf');


%---------------------------------------
load J:\matlab_hudson\probabilite_OUTPUT\annuelle\'Prob_2010.mat';
prob_8int=uint8(mean_prob_year);

cd('J:\matlab_hudson\probabilite_OUTPUT\HDF\yearly\');   
imwrite(prob_8int,'Prob_2010.hdf');

