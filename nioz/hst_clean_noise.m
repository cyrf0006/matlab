function hst_clean_noise(hdfFile, sensorNo)

% function hst_hdfQuickLook(hdfFile, nSkip)
%
% usage ex:
%  hst_hdfQuickLook('myFilename.h5', 100)
% 


T = hdf5read(hdfFile, '/data/temperature');
T = T';
n = hdf5read(hdfFile, '/dims/time');
z = hdf5read(hdfFile, '/dims/depth');

I = find(T<-10);
T(I) = NaN;

nPts = 1500;
allIndex = 1:length(n);
nLines = floor(length(allIndex)/nPts);
I = find(allIndex>nPts*nLines);
allIndex(I) = [];

Ind = reshape(allIndex, nPts, nLines)';
T0 = T;
for i = 1:length(sensorNo)
    TVec = T(sensorNo(i),:);
    for j = 1:size(Ind,1)
        keyboard
        figure(1)
        clf
        hist(TVec(Ind(j, :)), 20);
        [N, X] = hist(TVec(Ind(j, :)), 20);
        [Y,I] = max(N);
        N(1:I) = Y;
        [Y,I] = min(N);
               
        figure(2)
        clf
        plot(TVec(Ind(j,:)), 'b'); 
        if I == 1
            disp('WARNING, Check this')
            keyboard
        else
           Iremove = find(TVec(Ind(j,:))>=X(I)); 
           TVec(Ind(j,Iremove)) = NaN;
        end
       
       hold on
       plot(TVec(Ind(j,:)), 'r');        
       hold off
       pause
    end
    
    T(sensorNo(i), :) = TVec;
end
