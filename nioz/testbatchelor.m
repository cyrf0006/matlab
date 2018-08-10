clear

load Testdata.txt

n = Testdata(:,1);
temp = Testdata(:,2);
clear Testdata;

timebin = 600/86400; %sec
timeVec = n(1)+timebin/2:timebin:n(end);

fs = 1; %Hz
u = ones(length(n),1)*0.01; %m/s


visc = 1e-6;
% $$$ ax = ang2acc(pitch); ay = ang2acc(roll);
% $$$ DENS = sw_dens(SBS, SBT, P);
% $$$ [nu, mu] = viscosity(SBS, SBT, DENS);
% $$$ 
% $$$ fs = round(fs);
% $$$ FS = round(FS);

%dtdx = diff(temp(I),1,1)./diff(u(I)*timebin,1,1); 

for i = 1:length(timeVec)
    I = find(n>=timeVec(i)-timebin/2 & n< timeVec(i)+timebin/2);

    
    dtdx = gradient(temp(I))./(u(I).*timebin); 

    [ps0, f0] = pwelch(dtdx, [], [], [], 1); % Raw Power Spectrum
                                          
    % J = find(f0 ~= 0); % remove zeros
    
    k0 = f0/nanmean(u(I));
    
    % ----------- chi calculation ------------------- %
    nois=noise(k0); %get the noise at each value of k      ---> FC: noise must be adjusted
    dk=k0(3)-k0(2);
    Dt=.00000014; % Heat molec. diffus.
    chi_i = 6*Dt*sum((ps0'-nois).*dk); % (Eqn. 9)

    figure(2)
    clf
    loglog(k0, ps0, 'k')
    hold on
    plot([10 500], [nois(1) nois(1)], 'r')
    title(sprintf('chi = %d', chi_i))

    % This is main likelyhood calculation
    [chi_i, eps_i, likelihoodratio,likelihood,kb,f11] =  fit_kb(k0,ps0, 0);


    % BAtchelor Spectrum
     y = batchSpectrum(k0, chi_i, kb);     
     loglog(k0, y, 'r')
    
     keyboard
     % check: 
     %http://journals.ametsoc.org/doi/abs/10.1175/2009JTECHO611.1
    
end

A = [ax ay az];
z = P;

shearFreqHz = fs;
pressFreqHz = FS;
zbin = 1;

% -> from quick_chi.m
zs=interp1([1:length(z)],z,[1:length(dtdz)]'*(pressFreqHz/shearFreqHz));
if isempty(zs)==1 % no shear data to process
   epsilon = [];
   zfft = [];
   return;
end

z1 = ceil(min(abs(z)));  % top of profile
z2 = floor(max(abs(z)));  % bottom of profile
zfft = z1+zbin/2:zbin:z2-zbin/2; %reg. z prof. for FFT

if isempty(zfft), % no shear data to process
   epsilon = [];
   zfft = [];
   return;
end % if

I=find(zs>=z1 & zs<=z2); % we keep only full bins
dtdz=dtdz(I);
zs = zs(I); 
A = A(I,:); % if shear is adjusted, size of A must be too



i = 83;
I = find(zs>zfft(i)-zbin/2 & zs<zfft(i)+zbin/2); % I = Imicrosctruc
Ifine = find(z>zfft(i)-zbin/2 & z<zfft(i)+zbin/2); % Ifinesctruc

if (isempty(I)==1)
    chi(i) = NaN;
    continue
end

II = find(isnan(dtdz(I))==1);

if length(II)>round(length(I)/4)
    disp('skip bin')
    chi(i) = NaN;
    continue % skip if more than 25% NaNs in the bin
else
    dtdz(I(II))=0; % pad NaNs with zeros
end
% -------------------------------------------------------- %

% --- Compute FFT & clean spectrum ---- %
[ps0, f0] = pwelch(dtdz(I), [], [], [], 512); % Raw Power Spectrum
J = find(f0 ~= 0); % remove zeros
ps0 = ps0(J);
f0 = f0(J);
% ------------------------------------- %

% ---- Conversion Freq -> Wavenumber ---- %
Wm = mean(W(Ifine)); % mean falling speed for this bin 
k0 = f0/Wm;
%---------------------------------------- %

keyboard

% ----------- chi calculation ------------------- %
nois=noise(k0); %get the noise at each value of k      ---> FC: noise must be adjusted
dk=k0(3)-k0(2);
Dt=.00000014; % Heat molec. diffus.
chi_i = 6*Dt*sum((ps0'-nois).*dk); % (Eqn. 9)

loglog(k0, ps0, 'k')
hold on
plot([10 500], [nois(1) nois(1)], 'r')
title(sprintf('chi = %d', chi_i))


% BAtchelor Spectrum
load(['eps_' fname_in])
kb = (eps1(i)./nanmean(nu(Ifine))./Dt^2).^.25;


i = 3
I = find(zs>zfft(i)-zbin/2 & zs<zfft(i)+zbin/2); % I = Imicrosctruc
Ifine = find(z>zfft(i)-zbin/2 & z<zfft(i)+zbin/2); % Ifinesctruc

if (isempty(I)==1)
    chi(i) = NaN;
    continue
end

II = find(isnan(dtdz(I))==1);

if length(II)>round(length(I)/4)
    disp('skip bin')
    chi(i) = NaN;
    continue % skip if more than 25% NaNs in the bin
else
    dtdz(I(II))=0; % pad NaNs with zeros
end
% -------------------------------------------------------- %

% --- Compute FFT & clean spectrum ---- %
[ps0, f0] = pwelch(dtdz(I), [], [], [], 512); % Raw Power Spectrum
J = find(f0 ~= 0); % remove zeros
ps0 = ps0(J);
f0 = f0(J);
% ------------------------------------- %

% ---- Conversion Freq -> Wavenumber ---- %
Wm = mean(W(Ifine)); % mean falling speed for this bin 
k0 = f0/Wm;
%---------------------------------------- %

% ----------- chi calculation ------------------- %
nois=noise(k0); %get the noise at each value of k      ---> FC: noise must be adjusted
dk=k0(3)-k0(2);
Dt=.00000014; % Heat molec. diffus.
chi_i = 6*Dt*sum((ps0'-nois).*dk); % (Eqn. 9)

figure(2)
clf
loglog(k0, ps0)
hold on
plot([10 500], [nois(1) nois(1)], 'r')
title(sprintf('chi = %d', chi_i))
