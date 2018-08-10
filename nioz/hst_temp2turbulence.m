function hst_temp2turbulence(tempFile, TSrelFile)


% usage ex: hst_temp2turbulence('temp_wholeSerie_skip1000.mat', 'TSrel_cubic.mat')
%           hst_temp2turbulence('temp_wholeSerie_skip10.mat', 'TSrel_cubic.mat')
%           hst_temp2turbulence('temp_wholeSerie_skip10.mat', 'temp_rho_relFile.mat')
%           hst_temp2turbulence('temp_wholeSerie_skip1000.mat', 'temp_rho_relFile.mat')
a = 5; 
%p = [.0029 -.0701 .6118 33.3481];

load(tempFile) 

T = Titp;
% check orientation
[zVec, I] = sort(zVec);
T = T(I,:);

%S = hst_temp2sal(T, TSrelFile); 
rho = hst_temp2sal(T, TSrelFile); 
[rhoSort, Isort] = sort(rho, 1);
if nanmean(nanmean(rho))<500;
    disp('T-rho rel: provide sigma-T?')
    disp('  -> use rho = rho+1000')
    rho = rho+1000;
end

N2 = nan(size(rho));
for i = 1:size(rho,2) % explicit calculation (by-pass preceeding)
    N2(:,i) = 9.81./nanmean(rho(:,i)).*gradient(rhoSort(:,i))./gradient(zVec);
    %    N2(:,i) = 9.81./rho(:,i).*gradient(rhoSort(:,i))./gradient(zVec);
end


%% Thorpe scale method
%dz=gradient(zVec);
%[rhoSort,I]=sort(rho,1);
%zIndex = [1:size(rho,1)]';

ThScale = nan(size(rho));
for i=1:size(rhoSort,2)
    ThScale(:,i) = zVec(Isort(:,i)) - zVec;
end
eps=0.64*ThScale.^2.*sqrt(N2).^3;
kz=0.128*ThScale.^2.*sqrt(N2);

figure(3)
clf
semilogy(timeVec, nanmean(eps,1), 'b')
hold on

load browseData % calculated in broweDataRockall...
timeVec = tijd+datenum(2012, 1, 1);
semilogy(timeVec, nanmean(eps,1), 'r')
keyboard



% Thorpe scale
d = nan(size(N2));
for i = 1:size(T,2)
    dd = abs(zVec - zVec(flipud(I(:,i))));
    d(:,i) = interp1(zVec, dd, p_mid(:,1));    
end

epsilon = .64.*d.^2.*sqrt(N2).^3;

% Load REORANGE output
OVT = load('reorange_ovt');

% --- Hans' stuff --- %
TT = Titp;
dz = zVec;
T0 = TT;

[mT,nT]=size(TT);
sa=35.3*ones(mT,nT);
% UNcomment with Matlab!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dplct10=gsw_CT_from_t(sa,TT,dz*1.01-10.1325);
TT=dplct10;
nf1=1;
nf2=nT;

drdT=0.157*1e-2;
T=T0;
% T=flipdim(T0,1);
[mt,nt]=size(T);
nf1=1;
nf2=nt;

deltaz=1.0;
[Tsort,I]=sort(T,1);
[tmp,J]=sort(I,1);
[time,Z]=meshgrid([1:size(T,2)],[1:size(T,1)]);
%   ThScale=(I-Z)*deltaz;
[dTdt,dTTdz]=gradient(Tsort,1,deltaz);
ddz=-gradient(dz);
for ii=1:nt
    dTdz(:,ii)=dTTdz(:,ii)./ddz;
    ThScale(:,ii)=(I(:,ii)-Z(:,ii)).*ddz;
end
kz=0.128*ThScale.^2.*(drdT*dTdz).^0.5;
flux=kz.*dTdz;
Nbv=(drdT*dTdz).^0.5;
nu=1.e-6;
eps=0.64*ThScale.^2.*Nbv.^3;
kzhi=1.6*abs(ThScale).*(nu*Nbv).^0.5;
epsthres=100*nu*Nbv.^2;

% Plot
dt = 1/24;
timeVecCst = round(min(timeVec)):dt:round(max(timeVec));
epsAve = nan(size(timeVecCst));
epsAveOvt = nan(size(timeVecCst));
epsAveHans = nan(size(timeVecCst));
for i = 1:length(timeVecCst)
    I = find(timeVec>timeVecCst(i)-dt & timeVec<=timeVecCst(i)+dt);
    E = epsilon(:,I);
    epsAve(i) = nanmean(E(:));

    E = eps(:,I);
    epsAveHans(i) = nanmean(E(:));
 


    I = find(OVT.OVT(:,1)>timeVecCst(i)-dt & OVT.OVT(:,1)<=timeVecCst(i)+dt);
    epsVec = a.*OVT.OVT(I,5).*OVT.OVT(I,8);
    epsAveOvt(i) = nanmean(epsVec);
end


figure(1)
clf
plot(timeVecCst, epsAve, 'k', 'linewidth', 2)
set(gca, 'yscale', 'log')
datetick('x', 1)
set(gca, 'xgrid', 'on')
xlabel('time')
ylabel('\epsilon (W kg^{-1})')

hold on
plot(timeVecCst, epsAveOvt, 'r', 'linewidth', 2)
plot(timeVecCst, epsAveHans, 'm', 'linewidth', 2)
legend('Fred calc.', 'Hans calc.', 'APEF calc.')

% Spectral analysis
freq = 1/(timeVecCst(2)-timeVecCst(1))/86400; % Hz
%U = detrend(U10m);

nx = max(size(epsAve)); % Hanning
na = 1;
w = hanning(floor(nx/na));

I = find(isnan(epsAveHans));
epsAveHans(I) = 0;
[ps, f] = pwelch(epsAveHans, w, 0, [], freq*86400); 

figure(2)
clf
loglog(f, ps, 'k', 'linewidth', 2)  
hold on
ylabel('PSD')
xlim([2e-1 1.5e1])
ylim([1e-16 3e-13])
plot_harmonics5

text(1/12.42*24*1, 2e-13, 'M_2', 'fontsize', 9, 'fontweight', 'bold', 'horizontalAlignment', 'center')
text(1/12.42*24*1.5, 2e-13, 'M_3', 'fontsize', 9, 'fontweight', 'bold', 'horizontalAlignment', 'center')
text(1/12.42*24*2, 2e-13, 'M_4', 'fontsize', 9, 'fontweight', 'bold', 'horizontalAlignment', 'center')
%text(1/12.42*24*3, 2e-13, 'M_6', 'fontsize', 9, 'fontweight', 'bold', 'horizontalAlignment', 'center')
%text(1/12.42*24*4, 2e-13, 'M_8', 'fontsize', 9, 'fontweight', 'bold', 'horizontalAlignment', 'center')
text(1/15.981*24, 2e-13,'f', 'fontsize', 9, 'fontweight', 'bold', 'horizontalAlignment', 'center')
text(1/23.93*24, 2e-13,'K_1', 'fontsize', 9, 'fontweight', 'bold', 'horizontalAlignment', 'center')


keyboard


