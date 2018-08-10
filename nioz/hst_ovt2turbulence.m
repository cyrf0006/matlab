function hst_ovt2turbulence(tempfile, TSrelFile, ovtFiles)
    
% function hst_ovt2turbulence(tempfile, ovtFiles)
%
% (to be run in /home/cyrf0006/research/NIOZ/RockallBank (for example)
%
% Inspired from APEF_M2.m
%
% usage ex: 
% hst_ovt2turbulence('temp_wholeSerie_skip1000', 'TSrel_cubic.mat', 'ovt_skip1000.list')
% NO-hst_ovt2turbulence('temp_wholeSerie_skip10', 'TSrel_cubic.mat', 'ovt_skip10.list')
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
[Tsort, I] = sort(T);
S = hst_temp2sal(Tsort, TSrelFile); 
SA = gsw_SA_from_SP(S, zVec, -15.1, 55.5);
CT = gsw_CT_from_t(SA,Tsort,zVec);
rho = gsw_rho(SA,CT,zVec);
%[rho_sort, I] = sort(rho,1);
%[N2, p_mid] = gsw_Nsquared(SA(I),CT(I),flipud(zVec));
[N2, p_mid] = gsw_Nsquared(SA,CT,zVec);

% $$$ % For Octave
% $$$ rho = sw_dens(S,T,zVec);
% $$$ [N2,vort,p_mid] = sw_bfrq(S,T,zVec);
% $$$ 

% overturns
fid = fopen(ovtFiles);
C = textscan(fid, '%s', 'delimiter', '\n');
list_ovt = char(C{1});

no_profiles = size(list_ovt, 1);
if no_profiles ~= length(timeVec);
    disp(['Problem, time vector not the same size as number of *.ovt ' ...
          'files']);
    return
end

% initialize OVT matrix
OVT = [];

count=1;
for i = 1:no_profiles

    disp(sprintf('profile %d / %d', i, no_profiles));

    % -------- Work on overturns ------ %
    fname_ovt = list_ovt(i, :); % remove white space
    I = find(fname_ovt==' ');
    fname_ovt(I)=[];

    % get overturns for this profile
    command = sprintf('cut -d" " -f2 %s > /tmp/tmp', fname_ovt);
    system(command);
    ind1 = load('/tmp/tmp');

    command = sprintf('cut -d" " -f3 %s > /tmp/tmp', fname_ovt);
    system(command);
    ind2 = load('/tmp/tmp');
    
    command = sprintf('cut -d" " -f4 %s > /tmp/tmp', fname_ovt);
    system(command);
    z1 = load('/tmp/tmp');
    
    command = sprintf('cut -d" " -f5 %s > /tmp/tmp', fname_ovt);
    system(command);
    z2 = load('/tmp/tmp');
    
    command = sprintf('cut -d" " -f6 %s > /tmp/tmp', fname_ovt);
    system(command);
    L_t = load('/tmp/tmp'); 

% $$$     command = sprintf('cut -d" " -f7 %s > /tmp/tmp', fname_ovt);
% $$$     system(command);
% $$$     ThFluct = load('/tmp/tmp'); % <kg/m3>, cf. rho'
% $$$     
% $$$     command = sprintf('cut -d" " -f9 %s > /tmp/tmp', fname_ovt);
% $$$     system(command);
% $$$     L_t_on_H = load('/tmp/tmp');
% $$$     H = L_t./L_t_on_H;
    
    command = sprintf('cut -d" " -f11 %s > /tmp/tmp', fname_ovt);
    system(command);
    APEF = load('/tmp/tmp');
    
    if isempty(L_t)
        continue
    end
      
    no_ovt = length(L_t);
    for j = 1:no_ovt
        % Few calculations (for Ra & Jb)
        Hz = z2(j) - z1(j);
        
        I = find(p_mid(:,i)>=z1(j) &  p_mid(:,i)<=z2(j));
        N = nanmean(sqrt(N2(I,i)));
                
        I = find(zVec>=z1(j) &  zVec<=z2(j));
        rho_sort = sort(rho(I,i));
        rho0 = nanmean(rho_sort);
        rho_p = rho_sort(end)-rho_sort(1); % density anomaly
        Jb = APEF(j).*N;
        Ra = (g*rho_p*Hz.^3)./(rho0*Kt*nu);
        Ts = (Hz.^2./Jb).^(1/3);
       
        %OVT contains: [time Z1 Z2, L_t APEF, L_t ovtSize, N(brunt-V), Jb(buoFlux), Raleigh]
        if isempty(OVT)==1
            %zbottom come from 'N2_matrix', thus zbottom(i)-z1(j) is hab
            OVT = [timeVec(i) z1(j) z2(j) L_t(j) APEF(j), Hz, N, Jb, Ra];   
        else
            OVT = [OVT; [timeVec(i) z1(j) z2(j) L_t(j) APEF(j), Hz, N, Jb, Ra]]; 
        end
    end                
    % ------------------------------- % 
end


save reorange_ovt.mat OVT zVec timeVec

disp('reorange_ovt.mat saved!')
disp(' -> type ''return'' to see a bit of what turbulence looks like')
keyboard



z1 = OVT(:,2);
z2 = OVT(:,3);

APEF = [];
Jb = [];
Ra = [];
Y = [];
X = [];
for i = 1:length(z1)
    z = [round(z1(i)):round(z2(i))]';
    ti = z*0+OVT(i,1);
    Jbi = z*0+OVT(i,8);
    APEFi = z*0+OVT(i,5);%/(z2(i)-z1(i));
    Rai = z*0+OVT(i,9);
    
    Y = [Y; z];    
    X = [X; ti];
    Jb = [Jb; Jbi];
    APEF = [APEF; APEFi];
    Ra = [Ra; Rai];
end



figure(1) % Buoyancy fluxes
clf
[XI, YI] = meshgrid(timeVec, zVec);
ZI = griddata(X, Y, log10(Jb), XI, YI);
I = find(sum(~isnan(Traw),2)==0);
ZI(I,:) = NaN; 
pcolor(XI, YI, ZI); shading interp
hold on
contour(timeVec, p_mid(:,1), N2,20, 'k') 
%contour(timeVec, zVec, T, 5, 'k')
hold off
caxis([-9 -5])
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
hold on
contour(timeVec, p_mid, N2, 50, 'k')
hold off
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
hold on
contour(timeVec, p_mid, N2, 50, 'k')
hold off
%contourf(XI, YI, ZI, 100, 'linestyle', 'none');
datetick
set(gca, 'ydir', 'reverse')
xlim([min(timeVec) max(timeVec)])
colorbar
ylabel('Depth (m)')
title('R_a')

% $$$ figure(4) % Buoyancy Freq.
% $$$ clf
% $$$ pcolor(); shading interp
% $$$ %contourf(XI, YI, ZI, 100, 'linestyle', 'none');
% $$$ datetick
% $$$ set(gca, 'ydir', 'reverse')
% $$$ xlim([min(timeVec) max(timeVec)])
% $$$ colorbar
% $$$ ylabel('Depth (m)')
% $$$ title('R_a')
