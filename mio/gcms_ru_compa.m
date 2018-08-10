clear all
close all

all_gcms = ...
[0.4928, 3.2; ... %sauma01
0.5251,	43.3; ... 
0.4992,	11.4; ...
0.4992,	8.7; ...
0.4898,	0.7; ...
0.0692,	0.7; ... %sauma02
0.0721,	0.8; ...
0.0997,	2.1; ...
0.1523,	6.5; ...
0.1694,	15.4; ...
0.0562,	12.8; ... %troll01
0.0552,	1.1; ...
0.0554,	13.8; ...
0.0533,	0.8; ...
0.3175,	250.8; ... %cedre01
1.2543,	1917.9; ...
0.7598,	1079.1; ...
0.5277,	734.8; ...
0.4241,	600.5; ...
0.3853,	384.4];


% Remove sauma01 (wrong prism, blank too high)
%all_gcms(1:5, :) = [];



%p = polyfit(all_gcms(6:end,2),all_gcms(6:end,1),1) 
% $$$ xfit = 1:2000;
% $$$ yfit = p(2) + xfit.*p(1);

% $$$ %scatter(all_gcms(1:5,2), all_gcms(1:5,1), '.k')
% $$$ %hold on
% $$$ plot(all_gcms(6:10,2), all_gcms(6:10,1), '.m', 'markerSize', 20)
% $$$ hold on
% $$$ plot(all_gcms(11:14,2), all_gcms(11:14,1), '.r', 'markerSize', 20)
% $$$ plot(all_gcms(15:end,2), all_gcms(15:end,1), '.b', 'markerSize', 20)
% $$$ plot(xfit, yfit, 'k', 'linewidth', 2)
% $$$ legend('sauma02', 'troll01', 'cedre01')
% $$$ xlabel('Phes (ng/L)')
% $$$ ylabel('Phe-like (RU)')



p = polyfit(all_gcms(6:end,1),all_gcms(6:end,2),1) 
xfit = 0:.1:1.3;
yfit = p(2) + xfit.*p(1);
p = polyfit(all_gcms(6:end,1),all_gcms(6:end,2),2) 
%yfit2 = p(3) + xfit.*p(2) + xfit.^2.*p(1) ;

figure(1)
plot(all_gcms(6:10,1), all_gcms(6:10,2), '.m', 'markerSize', 20)
hold on
plot(all_gcms(11:14,1), all_gcms(11:14,2), '.r', 'markerSize', 20)
plot(all_gcms(15:end,1), all_gcms(15:end,2), '.b', 'markerSize', 20)
plot(xfit, yfit, 'k', 'linewidth', 2)
%plot(xfit, yfit2, 'g', 'linewidth', 2)
legend('sauma02', 'troll01', 'cedre01')
ylabel('Phes gcms (ng/L)')
xlabel('Phe-like (RU)')

figure(2)
subplot(131)
plot(all_gcms(6:10,1), all_gcms(6:10,2), '.r', 'markerSize', 20)
p = polyfit(all_gcms(6:10,1),all_gcms(6:10,2),1); 
xfit = 0:.01:.2;
yfit = p(2) + xfit.*p(1);
hold on
plot(xfit, yfit, 'k', 'linewidth', 2)
xlim([0, .2])
ylim([0, 20])
title(sprintf('sauma02 / %4.4fx + %f', p(1), p(2)), 'fontWeight', 'bold')
ylabel('Phes gcms (ng/L)', 'fontWeight', 'bold')
xlabel('Phe-like (RU)', 'fontWeight', 'bold')

subplot(132)
plot(all_gcms(11:14,1), all_gcms(11:14,2), '.r', 'markerSize', 20)
p = polyfit(all_gcms(11:14,1),all_gcms(11:14,2),1);
xfit = 0:.01:.2;
yfit = p(2) + xfit.*p(1);
hold on
plot(xfit, yfit, 'k', 'linewidth', 2)
xlim([0, .2])
ylim([0, 20])
title(sprintf('troll01 / %4.4fx + %f', p(1), p(2)), 'fontWeight', 'bold')
ylabel('Phes gcms (ng/L)', 'fontWeight', 'bold')
xlabel('Phe-like (RU)', 'fontWeight', 'bold')


subplot(133)
plot(all_gcms(15:end,1), all_gcms(15:end,2), '.r', 'markerSize', 20)
p = polyfit(all_gcms(15:end,1),all_gcms(15:end,2),1);
xfit = 0:.1:1.2;
yfit = p(2) + xfit.*p(1);
hold on
plot(xfit, yfit, 'k', 'linewidth', 2)
xlim([0, 1.5])
ylim([0, 2000])
title(sprintf('cedre01 / %4.4fx + %f', p(1), p(2)), 'fontWeight', 'bold')
ylabel('Phes gcms (ng/L)', 'fontWeight', 'bold')
xlabel('Phe-like (RU)', 'fontWeight', 'bold')