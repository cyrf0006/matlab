function [timeVec, KEVec, KEVec_top, KEVec_bot] = hst_KEtimeserie(tempFile, TSrelFile, adcpFile)
    
% usage ex: [timeVec, KE, KE1, KE2] = hst_KEtimeserie('Rockall_2minAve.mat', 'temp_rho_relFile_3rd.mat', 'roc12.mat')
    

%% Few parameters
g = 9.81; % m/s2    
    

%% From T to density
load(tempFile);
load(TSrelFile);
rho = polyval(p, Tbin)+1000;
[rhoSort, I] = sort(rho,1);


%% Velocity to sensor resolution
vel =  velocity2hst(adcpFile, zVec, timeVec);


%% KE calculation
H = abs(zVec(end)-zVec(1));
dz = gradient(zVec,1,1);
KE = nan(size(timeVec));
KE1 = nan(size(timeVec));
KE2 = nan(size(timeVec));
[Y, I] = min(abs(zVec-(max(zVec)-H/2)));

for i = 1:length(timeVec);
    KE(i) = nansum(rhoSort(:,i).*vel(:,i).^2.*dz)/H;
    KE1(i) = nansum(rhoSort(1:I,i).*vel(1:I,i).^2.*dz(1:I))/(H/2); %top
    KE2(i) = nansum(rhoSort(I:end,i).*vel(I:end,i).^2.*dz(I:end))/(H/2);  % bottom
end
KEVec = KE;
KEVec_top = KE1;
KEVec_bot = KE2;


figure(1)
clf
plot(timeVec, KE, 'linewidth', 2)
hold on
plot(timeVec, KE1, 'r')
plot(timeVec, KE2, 'm')
legend('KE', 'KE_{top}', 'KE_{bot}')
datetick('x', 7)
YLIM = get(gca, 'ylim');
xlim([min(timeVec) max(timeVec)])
