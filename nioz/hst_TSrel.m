function hst_TSrel(ctdList, timeFile, TLIM)
    
% function hst_TSrel(ctdList, timeFile)
%
% Made specificly to read and process CTD casts from ROC12
% cruise. Shell preccessing were done to clean file names and to
% build the timeFile
%    
% usage ex: hst_TSrel('ctdFiles.list','./CTD360/fred_processed/timeFile.txt', [])



% load CTD files
fid = fopen(ctdList);
C = textscan(fid, '%s', 'delimiter', '\n');
ctd = char(C{1});
noFiles = size(ctd,1);

% get time vector
fid = fopen(timeFile);
C = textscan(fid, '%s', 'delimiter', '\n');
ctdTime = char(C{1});

TVec = [];
ZVec = [];
SVec = [];
timeVec = [];

for i = 1:noFiles

    fname = ctd(i, :);
    I = find(fname==' ');   
    fname(I) = [];
    
    data = load(fname);
   
    Z = data(:, 14);
    T = data(:, 15);
    S = data(:, 16);
    
     % Isolate downcast
    dz = gradient(Z);
    I = find(dz>0.5*mean(abs(dz)));
    
    ZVec = [ZVec; Z(I)];
    TVec = [TVec; T(I)];
    SVec = [SVec; S(I)];
    timeVec = [timeVec; Z(I)*0+ctdTime(i)];
    
end

% Clean T,S
I = find(SVec<1e-6);
SVec(I) = [];
TVec(I) = [];
ZVec(I) = [];
I = find(TVec<1e-6);
SVec(I) = [];
TVec(I) = [];
ZVec(I) = [];


figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 10])
I = find(TVec>=TLIM(1) & TVec <= TLIM(2));
TVec = TVec(I);
SVec = SVec(I);


[Y, I] = sort(TVec);
p = polyfit(TVec(I), SVec(I), 1);
Tfit = min(TVec):.1:max(TVec);
Sfit = Tfit*p(1) + p(2);


p = polyfit(TVec(I), SVec(I), 3);
Sfit3 = Tfit.^3*p(1) + Tfit.^2*p(2) + Tfit*p(3) + p(4);


plot(TVec, SVec, '.k')
hold on
plot(Tfit, Sfit, 'r') 
plot(Tfit, Sfit3, 'm') 


% $$$ keyboard
% $$$  Sfit = Y*p(1) + p(2);
% $$$ Sfit3 = Y.^3*p(1) + Y.^2*p(2) + Y*p(3) + p(4);
% $$$ corrcoef(Tfit, Y)
% $$$ corrcoef(Tfit3, Y)
