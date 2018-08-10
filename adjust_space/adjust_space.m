% adjust_space.m
% 
% This scripts is supose to automaticly adjust space between
% subplots to minimize blank or wathever between them. To us it,
% just copy and edit the following lines at the beginning of any
% script that plot subplots. Then, after you've call and custom
% your plot, just call "adjust_space". The script should do the
% rest!
%
% Lines to copy at the beginning of any scripts:
%

% $$$ % *********************** Adjust_space.m ************************ %
% $$$ % Fields required by the function adjust_space.m. Please fill every
% $$$ % of the following and call "adjust_space" in the script whenever
% $$$ % you want. Do not touch four last fields
% $$$ ncol = 4; % no. subplot column
% $$$ nrow = 2; % no. subplot row
% $$$ dx = 0.03 ; % horiz. space between subplots
% $$$ dy = 0.04; % vert. space between subplots
% $$$ lefs = 0.1; % very left of figure
% $$$ rigs = 0.1; % very right of figure
% $$$ tops = 0.05; % top of figure
% $$$ bots = 0.1; % bottom of figure
% $$$ figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
% $$$ figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
% $$$ count_col = 1;
% $$$ count_row = 1;
% $$$ % *************************************************************** %
%
% NOTE: When there is more than one figure to adjust, please call
% adjust_space.m only for one at a time. This could be easy to
% correct, because it comes from count_row and count_col. You must
% not increment them if there is another figure to adjust!
%
% Author: Frederic Cyr - April 2011
%
% --------------------------------------------------------------- %



pos = get(gca, 'position'); % only to fill it

pos(3) = figw; % figure width
pos(4) = figh; % figure heigth
pos(1) = lefs+(count_col-1)*(figw+dx); % x0
pos(2) = 1-tops-count_row*figh-(count_row-1)*dy; % y0

count_col = count_col + 1;

if count_col > ncol
    count_col = 1;
    count_row = count_row + 1;
end

set(gca, 'pos', pos)

