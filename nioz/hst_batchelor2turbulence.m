function hst_batchelor2turbulence(hdfFile, adcpFile)

% function hst_hdf2turbulence(hdfFile)
%
% usage ex:
%  hst_batchelor2turbulence('myFilename.h5', '/media/Seagate1TB/NIOZ/thermistordata/ROC12/roc12.mat')
% 
% legend('tbin600-fft100', 'tbin1200-fft100', 'tbin1200-fft300', 'tbin300-fft30')

%% Few parameters
visc = 1e-6; 
NIOZ3_noise = 1e-4; % 0.1mK
nois = 0;
fs = 1; %Hz
tbin = 300; %sec
%tbin = 300; %sec
fftbin = 30; % sec

% $$$ tbin = 10; 
% $$$ fftbin = 10;

timebin = tbin/86400; %in day

%% read HDF file
T = hdf5read(hdfFile, '/data/temperature');
T = T';
n = hdf5read(hdfFile, '/dims/time'); % raw sensor time
zVec = hdf5read(hdfFile, '/dims/depth');

% new timeVector
timeSensor = n+datenum(2012, 0, 0); % THIS MUST BE GENERALIZED!!!!!!!!(EX:2012)
clear n
t0 = timeSensor(1); %n(1)+datenum(2012, 0, 0);
tf = timeSensor(end); %n(end)+datenum(2012, 0, 0);
%timeVec = n(1)+timebin/2:timebin:n(end);
timeVec = t0+timebin/2:timebin:tf;


%% Read ADCP file
% assuming constant velocity for now
%u = ones(1, length(n))*0.01; %m/s
uAdcp = velocity2hst(adcpFile, zVec, timeVec); 


%% Process data
epsMat = nan(length(zVec), length(timeVec));
chiMat = nan(length(zVec), length(timeVec));
kbMat = nan(length(zVec), length(timeVec));

for iz = 1:length(zVec)
    disp(sprintf('sensor %d / %d', iz, length(zVec)))
    temp = T(iz,:);

    if sum(~isnan(temp)) == 0
        disp(' -> bad sensor (continue)')
        continue
    end
    
    for i = 1:length(timeVec)
        I = find(timeSensor>=timeVec(i)-timebin/2 & timeSensor< timeVec(i)+timebin/2);        
        
% $$$         if mod(i,25) == 0
% $$$             disp(sprintf('%d', i))
% $$$         end
% $$$         
        
        if isnan(uAdcp(iz,i))
            continue
        end
    
        dtdx = gradient(temp(I))./(uAdcp(iz,i).*tbin);    
        II = find(isnan(dtdx));
        if ~isempty(II)
            dtdx(II) = 0;
        end
        
        % if detrend is needed
% $$$     ddtdx = detrend(dtdx);
% $$$     [ps1, f1] = pwelch(ddtdx, [], [], [], fs); % Raw Power Spectrum
% $$$     ps1 = ps1*nanmean(u(I));
% $$$     k1 = f1/nanmean(u(I));

        nx = max(size(dtdx)); % Hanning   
        na = round(tbin/fftbin);
        w = hanning(floor(nx/na));
        
        [ps0, f0] = pwelch(dtdx, w, [], [], fs); % Raw Power Spectrum    
        ps0 = ps0*nanmean(uAdcp(iz,i)); % (degC/m)^2/Hz -> (degC/m)^2/cpm 
        k0 = f0/nanmean(uAdcp(iz,i)); % Hz -> cpm
        
        % ----------- chi calculation ------------------- %
        %nois=noise(k0); %get the noise at each value of k      ---> FC: noise must be adjusted
        dk=k0(3)-k0(2);
        Dt=.00000014; % Heat molec. diffus.
        chi_i = 6*Dt*sum((ps0'-nois).*dk); % (Eqn. 9)
        
        % This is main likelyhood calculation
        [chi_i, eps_i, likelihoodratio, likelihood, kb, f11] =  fit_kb(k0,ps0, 0);

        % BAtchelor Spectrum
        y = batchSpectrum(k0, chi_i, kb);     
        
        
    figure(2)
        clf
        loglog(k0, ps0, 'k')
        hold on
        %plot([10 500], [nois(1) nois(1)], 'r')
        loglog(k0, y, 'r')
        title(sprintf('chi = %d; eps = %d; kb = %d', chi_i, eps_i, kb))
        keyboard
        
        % check: 
        %http://journals.ametsoc.org/doi/abs/10.1175/2009JTECHO611.1
        
        epsMat(iz, i) = eps_i;
        chiMat(iz, i) = chi_i;
        kbMat(iz, i) = kb;
        
        
    end
end

keyboard


% $$$ 
% $$$ i=1110; iz=30;
% $$$ i=923; iz = 29;
% $$$ i=2033; iz = 30;
% $$$ i=1114; iz=17;
% $$$ temp = T(iz,:);
% $$$ I = find(timeSensor>=timeVec(i)-timebin/2 & timeSensor< timeVec(i)+timebin/2);        
% $$$ dtdx = gradient(temp(I))./(uAdcp(iz,i).*tbin);    
% $$$ II = find(isnan(dtdx));
% $$$ if ~isempty(II)
% $$$     dtdx(II) = 0;
% $$$ end
% $$$ 
% $$$ nx = max(size(dtdx)); % Hanning   
% $$$ na = round(tbin/fftbin);
% $$$ w = hanning(floor(nx/na));
% $$$ 
% $$$ [ps0, f0] = pwelch(dtdx, w, [], [], fs); % Raw Power Spectrum    
% $$$ ps0 = ps0*nanmean(uAdcp(iz,i)); % (degC/m)^2/Hz -> (degC/m)^2/cpm 
% $$$ k0 = f0/nanmean(uAdcp(iz,i)); % Hz -> cpm
% $$$ 
% $$$ dk=k0(3)-k0(2);
% $$$ Dt=.00000014; % Heat molec. diffus.
% $$$ chi_i = 6*Dt*sum((ps0'-nois).*dk); % (Eqn. 9)
% $$$ [chi_i, eps_i, likelihoodratio, likelihood, kb, f11] =  fit_kb(k0,ps0, 0);
% $$$ y = batchSpectrum(k0, chi_i, kb);     
% $$$ 
% $$$ figure(2)
% $$$ clf
% $$$ loglog(k0, ps0, 'k')
% $$$ hold on
% $$$ loglog(k0, y, 'r')
% $$$ title(sprintf('chi = %d; eps = %d; kb = %d', chi_i, eps_i, kb))
