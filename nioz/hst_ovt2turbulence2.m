function hst_ovt2turbulence2(tempfile, TSrelFile)
    
% function hst_ovt2turbulence(tempfile, ovtFiles)
%
% (to be run in /home/cyrf0006/research/NIOZ/RockallBank (for example)
%
% Inspired from APEF_M2.m
%
% usage ex: 
% hst_ovt2turbulence2('temp_wholeSerie_skip1000', 'TSrel_cubic.mat')
% hst_ovt2turbulence2('temp_wholeSerie_skip1000', 'temp_rho_relFile.mat') <-- uses temp-rho relationship
%
% *** in a shell do (to generate *.ovt):
%   reorange -c 1 *.dat; 
%
% F. Cyr - April 2014
    
% few params:
g = 9.81; %m/s2
nu = 1e-6;% m2/s
Kt = 1.4e-7; %m2/s

% profiles
load(tempfile);

T = Titp;

% Compute N2 (from S)
% $$$ [Tsort, I] = sort(T);
% $$$ S = hst_temp2sal(Tsort, TSrelFile); 
% $$$ SA = gsw_SA_from_SP(S, zVec, -15.1, 55.5);
% $$$ CT = gsw_CT_from_t(SA,Tsort,zVec);
% $$$ rho = gsw_rho(SA,CT,zVec);

% Compute N2 (directly rho)
% $$$ [Tsort, I] = sort(T);
% $$$ rho = hst_temp2sal(Tsort, TSrelFile); 
% $$$ [N2, p_mid] = gsw_Nsquared(SA,CT,zVec);
% $$$ 
% $$$ N2 = buoy_freq(rho, zVec, zVec);
% $$$ N2 = -g./rho.*gradient(rho)./gradient(zVec);

% Re-compute T,S,rho (unsorted)
% $$$ S = hst_temp2sal(T, TSrelFile); 
% $$$ SA = gsw_SA_from_SP(S, zVec, -15.1, 55.5);
% $$$ CT = gsw_CT_from_t(SA,T,zVec);
% $$$ rho = gsw_rho(SA,CT,zVec);

rho = hst_temp2sal(T, TSrelFile); 

% make sure the profile is not upside down:
[zVec, I] = sort(zVec);
rho = rho(I,:);
T = T(I,:);


% initialize OVT matrix
OVT = [];

count=1;

for i = 1:size(T,2)

    disp(sprintf('profile %d / %d', i, size(T,2)));

    rho_raw = rho(:,i);    
    [rho_sort, I] = sort(rho_raw);
% $$$     plot(rho_raw, zVec)
% $$$     hold on
% $$$     plot(rho_sort, zVec, 'r')
% $$$     set(gca, 'ydir', 'reverse')

    % find overturning sequences
    J = find(I-[1:length(I)]'==0);
        
    % clean J
    %   J = [1; J; length(I)];
    II = find(diff(J)<3); %stable portion at least length 2
    J = J(II);
    
    J = [1; J; length(I)];
    
    limits = [];
    
    if isempty(J) % no ovt
        continue
    elseif length(J) == 1 % one ovt (whole profile)
        limits = [1 length(I)];
        theEnd = 1;
    else % at least 2 ovt
       theEnd = 0; 
       iovt=1;
    end        

    while ~theEnd                
        
        if J(iovt+1)-J(iovt)>3 % unstable portion must be at least length=3
            limits = [limits; J(iovt)+1 J(iovt+1)]; 
            iovt = iovt+1;
        else
            iovt = iovt+1;
        end
        
        if iovt == length(J);
            theEnd = 1;
        end
    end
    
    % correct if index(1) is two
    if limits(1,1) == 2
        limits(1,1) = 1;
    end
            
    % loop on ovt
    no_ovt = size(limits,1);

    
    APEF = [];
    Jb = [];
    Ra = [];
    Y = [];
    X = [];
    for j = 1:no_ovt
        % Few calculations (for Ra & Jb)

        I1 = limits(j,1);
        I2 = limits(j,2);
        z2 = zVec(I2);
        z1 = zVec(I1);
        
        Hz = zVec(I2)-zVec(I1); % abs cause zVec upsidedown
                
        rho_raw_ovt = rho_raw(I1:I2);
        [rho_sort_ovt, I]= sort(rho_raw_ovt);        
        rho0 = nanmean(rho_raw_ovt);
        zOvt = zVec(I1:I2);

        N = g/rho0*gradient(rho_sort_ovt)./gradient(zOvt);
        d = zOvt - zOvt(I);
        rho_p = rho_raw_ovt-rho_sort_ovt;
        
        ksi = g/rho0/Hz.*sum(rho_p.*(max(zOvt)-zOvt).*gradient(zOvt));
        
        %        APEF = [APEF; ksi];
                
        Jb = ksi.*nanmean(N);
        Ra = (g*(max(rho_p))*Hz.^3)./(rho0*Kt*nu); %double check rho_p
        Ts = (Hz.^2./Jb).^(1/3);

        e1 = .64*rms(d).^2*nanmean(N).^3;
% $$$         if e1<1e-15
% $$$             keyboard
% $$$         end
% $$$         e2 = nanmean(.64*d.^2*nanmean(N).^3);
% $$$         if abs(e1-e2) > 10
% $$$             disp('!')
% $$$             keyboard
% $$$         end

        %OVT contains: [time Z1 Z2, L_t APEF, L_t ovtSize, N(brunt-V), Jb(buoFlux), Raleigh]
        if isempty(OVT)==1
            %zbottom come from 'N2_matrix', thus zbottom(i)-z1(j) is hab
            OVT = [timeVec(i), z1, z2, rms(d), ksi, Hz, nanmean(N), ...
                   Jb, Ra, e1];   
        else
            OVT = [OVT; [timeVec(i), z1, z2, rms(d), ksi, Hz, ...
                         nanmean(N), Jb, Ra, e1]]; 
        end
    end                
    % ------------------------------- % 
end


save mooring_ovt.mat OVT zVec timeVec

disp('reorange_ovt.mat saved!')
disp(' -> type ''return'' to see a bit of what turbulence looks like')

z1 = OVT(:,2);
z2 = OVT(:,3);

APEF = [];
Jb = [];
Ra = [];
Y = [];
X = [];
EPS = [];

for i = 1:length(z1)
    z = [round(z1(i)):round(z2(i))]';
    ti = z*0+OVT(i,1);
    Jbi = z*0+OVT(i,8);
    APEFi = z*0+OVT(i,5);
    Rai = z*0+OVT(i,9);
    EPSi = z*0+OVT(i,10);

% $$$     N = OVT(i,7);
% $$$     RMSi = OVT(i,4);
% $$$     EPSi = z*0+.64*RMSi.^2.*N.^3;
        
    Y = [Y; z];    
    X = [X; ti];
    Jb = [Jb; Jbi];
    APEF = [APEF; APEFi];
    Ra = [Ra; Rai];
    EPS = [EPS; EPSi];
end

keyboard


figure(1) % Buoyancy fluxes
clf
[XI, YI] = meshgrid(timeVec, zVec);
ZI = griddata(X, Y, log10(Jb), XI, YI);
I = find(sum(~isnan(Traw),2)==0);
ZI(I,:) = NaN; 
pcolor(XI, YI, ZI); shading interp
cb = colorbar;
CAXIS = get(gca, 'clim');
CAXIS(1) = CAXIS(1)+2;
hold on
contour(timeVec, zVec, rho,20, 'k') 
%contour(timeVec, zVec, T, 5, 'k')
hold off
caxis([CAXIS])
%caxis([-9 -5])
%contourf(XI, YI, ZI, 100, 'linestyle', 'none');
datetick
xlim([min(timeVec) max(timeVec)])
set(gca, 'ydir', 'reverse')
colorbar
ylabel('Depth (m)')
title('J_b')

figure(2) % APEF
clf
[XI, YI] = meshgrid(timeVec, zVec);
ZI = griddata(X, Y, log10(APEF), XI, YI);
I = find(sum(~isnan(Traw),2)==0);
ZI(I,:) = NaN; 
pcolor(XI, YI, ZI); shading interp
CAXIS = get(gca, 'clim');
CAXIS(2) = CAXIS(2)+2;
hold on
contour(timeVec, zVec, rho,20, 'k') 
hold off
caxis([CAXIS])
%contourf(XI, YI, ZI, 100, 'linestyle', 'none');
datetick
xlim([min(timeVec) max(timeVec)])
set(gca, 'ydir', 'reverse')
colorbar
ylabel('Depth (m)')
title('\ksi')

figure(3) % Raleigh Number
clf
[XI, YI] = meshgrid(timeVec, zVec);
ZI = griddata(X, Y, log10(Ra), XI, YI);
I = find(sum(~isnan(Traw),2)==0);
ZI(I,:) = NaN; 
pcolor(XI, YI, ZI); shading interp
CAXIS = get(gca, 'clim');
hold on
CAXIS(2) = CAXIS(2)+2;
contour(timeVec, zVec, rho,20, 'k') 
hold off
caxis([CAXIS])
%contourf(XI, YI, ZI, 100, 'linestyle', 'none');
datetick
set(gca, 'ydir', 'reverse')
xlim([min(timeVec) max(timeVec)])
colorbar
ylabel('Depth (m)')
title('R_a')

figure(4) % Buoyancy Freq.
clf
pcolor(); shading interp
%contourf(XI, YI, ZI, 100, 'linestyle', 'none');
datetick
set(gca, 'ydir', 'reverse')
xlim([min(timeVec) max(timeVec)])
colorbar
ylabel('Depth (m)')
title('R_a')

keyboard

load browseData    
timeVec2 = tijd+datenum(2012, 1, 1);
epsThorpe = nanmean(eps,1);                               
[XI, YI] = meshgrid(timeVec, zVec);
ZI = griddata(X, Y, Jb./.2, XI, YI);
EPS = griddata(X, Y, EPS, XI, YI);
epsAPEF = nanmean(ZI,1); 
epsRMS = nanmean(EPS);
epsd = nanmean(eps,1);
figure(4)
clf
semilogy(timeVec, epsAPEF, 'b')
hold on
semilogy(timeVec, epsRMS, 'm')
semilogy(timeVec2, epsd, 'r')
ylabel('\epsilon (W kg^{-1})')
datetick
legend('APEF=1.5e-6','Thorpe = 3e-7')   



figure(5) % spectra
clf
freq1 = 1/(timeVec(2)-timeVec(1)); % pd
freq2 = 1/(timeVec2(2)-timeVec2(1)); 

nx = max(size(timeVec)); % Hanning
na = 1;
w = hanning(floor(nx/na));
I = find(isnan(epsAPEF));
epsAPEF(I) = 0;
I = find(isnan(epsRMS));
epsRMS(I) = 0;
[ps1, f1] = pwelch(epsAPEF, w, 0, [], freq1/86400); 
[ps2, f2] = pwelch(epsRMS, w, 0, [], freq1/86400); 

nx = max(size(timeVec2)); % Hanning
na = 1;
w = hanning(floor(nx/na));
I = find(isnan(epsd));;
epsd(I) = 0;
[ps3, f3] = pwelch(epsd, w, 0, [], freq2/86400); 

loglog(f1*86400, ps1)
hold on
loglog(f2*86400, ps2, 'm')
loglog(f3*86400, ps3, 'r')
ylabel('PSD')
xlabel('cpd')
xlim([1e-1 1e1])
%ylim([1e-22 3e-12])
plot_harmonics5
