% Example script: Copy-paste the following into a script to try it!
% ############## EXAMPLE SCRIPT BEGINNING ##################
figure(1)
clf
set(gcf,'PaperUnits','centimeters', 'paperposition', [0 0 20 20])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 2; % no. subplot column
nrow = 3; % no. subplot row
dx = 0.1 ; % horiz. space between subplots
dy = 0.1; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.1; % very right of figure
tops = 0.1; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

for i = 1:(ncol.*nrow)

    subplot(nrow,ncol,i)
    
    hist(rand(100,1))
    
    if i == 2 | i == 4 | i == 6
        set(gca, 'yticklabel', '')
    end
    if i == 5 
        xlabel('xlabel')
    end
    if i == 1 | i == 2 | i == 3 | i == 4
        set(gca, 'xticklabel', '')
    end
    if i == 3
        ylabel('Freq')
    end
    
    adjust_space
end
% ############# EXAMPLE SCRIPT END ####################
%