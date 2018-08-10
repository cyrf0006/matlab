clear

% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.0 ; % horiz. space between subplots
dy = 0.0; % vert. space between subplots
lefs = 0.05; % very left of figure
rigs = 0.05; % very right of figure
tops = 0.05; % top of figure
bots = 0.05; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

paperwidth = 16;%cm
paperheight = 16;%cm

figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 paperwidth paperheight])
%set(gcf, 'renderer', 'opengl')
m_proj('stereographic','latitude',60, 'longitude', -100, 'radius',52,'rotangle',0);
m_coast('patch','k');
m_grid('linewi',2,'tickdir','out');
%m_gshhs_l('patch',[.7 .7 .7]); %coastlines                               
%m_patch([0 .5333 .6667]) % blue gebco
m_gshhs_h('patch',[1 .9333 .6667]); %coastlines (Beige)                               


% Regions
%lon_vec = [-77:-35];
%lat_vec = [30:66];
% $$$ lon_vec = [-77:-51];
% $$$ lat_vec = [40:66];
% $$$ lon_atlantic = [lon_vec, fliplr(lon_vec), lon_vec(1)];
% $$$ lat_atlantic = [ones(size(lon_vec))*lat_vec(1), ones(size(lon_vec))*lat_vec(end), lat_vec(1)];


% $$$ lon_vec = [-79:-60];
% $$$ lat_vec = [58:65];
% $$$ lon_hudSt = [lon_vec, fliplr(lon_vec), lon_vec(1)];
% $$$ lat_hudSt = [ones(size(lon_vec))*lat_vec(1), ones(size(lon_vec))*lat_vec(end), lat_vec(1)];
% $$$ 
% $$$ lon_vec = [-63:-54];
% $$$ lat_vec = [52:59];
% $$$ lon_labShelf = [lon_vec, fliplr(lon_vec), lon_vec(1)];
% $$$ lat_labShelf = [ones(size(lon_vec))*lat_vec(1), ones(size(lon_vec))*lat_vec(end), lat_vec(1)];

lon_vec = [-79:-54];
lat_vec = [52:65];
lon_hudStLab = [lon_vec, fliplr(lon_vec), lon_vec(1)];
lat_hudStLab = [ones(size(lon_vec))*lat_vec(1), ones(size(lon_vec))*lat_vec(end), lat_vec(1)];

lon_vec = [-72:-58];
lat_vec = [41:45];
lon_nsShelf = [lon_vec, fliplr(lon_vec), lon_vec(1)];
lat_nsShelf = [ones(size(lon_vec))*lat_vec(1), ones(size(lon_vec))*lat_vec(end), lat_vec(1)];

lon_vec = [-60:-44];
lat_vec = [44:52];
lon_nfdl = [lon_vec, fliplr(lon_vec), lon_vec(1)];
lat_nfdl = [ones(size(lon_vec))*lat_vec(1), ones(size(lon_vec))*lat_vec(end), lat_vec(1)];

lon_vec = [-70:-55];
lat_vec = [45:52];
lon_gsl = [lon_vec, fliplr(lon_vec), lon_vec(1)];
lat_gsl = [ones(size(lon_vec))*lat_vec(1), ones(size(lon_vec))*lat_vec(end), lat_vec(1)];

lon_vec = [-85:-50];
lat_vec = [60:82];
lon_baffin = [lon_vec, fliplr(lon_vec), lon_vec(1)];
lat_baffin = [ones(size(lon_vec))*lat_vec(1), ones(size(lon_vec))*lat_vec(end), lat_vec(1)];

lon_vec = [-96:-72];
lat_vec = [51:71];
lon_hudson = [lon_vec, fliplr(lon_vec), lon_vec(1)];
lat_hudson = [ones(size(lon_vec))*lat_vec(1), ones(size(lon_vec))*lat_vec(end), lat_vec(1)];

lon_vec = [-180:-110];
lat_vec = [40:67];
lon_pacific = [lon_vec, fliplr(lon_vec), lon_vec(1)];
lat_pacific = [ones(size(lon_vec))*lat_vec(1), ones(size(lon_vec))*lat_vec(end), lat_vec(1)];

lon_vec = [-145:-115];
lat_vec = [67:76];
lon_beaufort = [lon_vec, fliplr(lon_vec), lon_vec(1)];
lat_beaufort = [ones(size(lon_vec))*lat_vec(1), ones(size(lon_vec))*lat_vec(end), lat_vec(1)];


% $$$ %m_patch(lon_atlantic, lat_atlantic , [1 .9333 .6667*1.05]*1.05, 'linewidth', 2)
% $$$ m_patch(lon_hudStLab, lat_hudStLab , [1 .9333*1.05 .6667*1.05], 'linewidth', 2)
% $$$ m_patch(lon_nsShelf, lat_nsShelf , [1 .9333*1.05 .6667*1.05], 'linewidth', 2)
% $$$ m_patch(lon_nfdl, lat_nfdl , [1 .9333*1.05 .6667*1.05], 'linewidth', 2)
% $$$ m_patch(lon_gsl, lat_gsl , [1 .9333*1.05 .6667*1.05], 'linewidth', 2)
% $$$ m_patch(lon_baffin, lat_baffin , [1 .9333*1.05 .6667*1.05], 'linewidth', 2)
% $$$ m_patch(lon_hudson, lat_hudson , [1 .9333*1.05 .6667*1.05], 'linewidth', 2)
% $$$ m_patch(lon_pacific, lat_pacific , [1 .9333*1.05 .6667*1.05], 'linewidth', 2)
% $$$ m_patch(lon_beaufort, lat_beaufort , [1 .9333*1.05 .6667*1.05], 'linewidth', 2, 'linestyle', '--')
% $$$ m_hatch(lon_beaufort, lat_beaufort,'single',30,5,'color','k'); % ...with hatching added.

m_patch(lon_hudStLab, lat_hudStLab , [1 .9333 .6667]*.9, 'linewidth', 2)
m_patch(lon_nsShelf, lat_nsShelf , [1 .9333 .6667]*.9, 'linewidth', 2)
m_patch(lon_nfdl, lat_nfdl , [1 .9333 .6667]*.9, 'linewidth', 2)
m_patch(lon_gsl, lat_gsl , [1 .9333 .6667]*.9, 'linewidth', 2)
m_patch(lon_baffin, lat_baffin , [1 .9333 .6667]*.9, 'linewidth', 2)
m_patch(lon_hudson, lat_hudson , [1 .9333 .6667]*.9, 'linewidth', 2)
m_patch(lon_pacific, lat_pacific , [1 .9333 .6667]*.9, 'linewidth', 2)
m_patch(lon_beaufort, lat_beaufort , [1 .9333 .6667]*.9, 'linewidth', 2, 'linestyle', '--')
m_hatch(lon_beaufort, lat_beaufort,'single',30,5,'color','k'); % ...with hatching added.



%m_text(max(lon_atlantic)-3, min(lat_atlantic)+2, '(A)', 'color', 'k','FontSize',6, 'fontweight', 'bold', 'BackgroundColor',[1 1 1], 'vertical', 'middle', 'horizontal', 'center', 'rotation', 45)

m_text(min(lon_hudStLab)+4, min(lat_hudStLab)+2.5, '(A)', 'color', 'k','FontSize',6, 'fontweight', 'bold', 'BackgroundColor',[1 1 1], 'vertical', 'middle', 'horizontal', 'center', 'rotation', 24.5)

m_text(min(lon_nfdl)+3, min(lat_nfdl)+2, '(B)', 'color', 'k','FontSize',6, 'fontweight', 'bold', 'BackgroundColor',[1 1 1], 'vertical', 'middle', 'horizontal', 'center', 'rotation', 35)

m_text(min(lon_nsShelf)+3, min(lat_nsShelf)+2, '(C)', 'color', 'k','FontSize',6, 'fontweight', 'bold', 'BackgroundColor',[1 1 1], 'vertical', 'middle', 'horizontal', 'center', 'rotation', 28.5)

m_text(min(lon_gsl)+3, min(lat_gsl)+2, '(D)', 'color', 'k','FontSize',6, 'fontweight', 'bold', 'BackgroundColor',[1 1 1], 'vertical', 'middle', 'horizontal', 'center', 'rotation', 30)

m_text(max(lon_baffin)-6, min(lat_baffin)+2, '(F)', 'color', 'k','FontSize',6, 'fontweight', 'bold', 'BackgroundColor',[1 1 1], 'vertical', 'middle', 'horizontal', 'center', 'rotation', 35)

m_text(min(lon_hudson)+4, min(lat_hudson)+2.5, '(E)', 'color', 'k','FontSize',6, 'fontweight', 'bold', 'BackgroundColor',[1 1 1], 'vertical', 'middle', 'horizontal', 'center', 'rotation', 5)

m_text(max(lon_pacific)-4, min(lat_pacific)+2.5, '(G)', 'color', 'k','FontSize',6, 'fontweight', 'bold', 'BackgroundColor',[1 1 1], 'vertical', 'middle', 'horizontal', 'center', 'rotation', -14)

m_text(min(lon_beaufort)+8, max(lat_beaufort)-2.5, '(H)', 'color', 'k','FontSize',6, 'fontweight', 'bold', 'BackgroundColor',[1 1 1], 'vertical', 'middle', 'horizontal', 'center', 'rotation', -37.6)

adjust_space


set(gcf, 'renderer', 'painters')
print('-dpng', '-r300', 'global_map.png')
%print('-depsc2', 'global_map.eps')
