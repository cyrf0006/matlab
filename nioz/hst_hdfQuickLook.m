function hst_hdfQuickLook(hdfFile, nSkip)

% function hst_hdfQuickLook(hdfFile, nSkip)
%
% usage ex:
%  hst_hdfQuickLook('myFilename.h5', 100)
% 


% read HDF file
T = hdf5read(hdfFile, '/data/temperature');
T = T';
n = hdf5read(hdfFile, '/dims/time');
z = hdf5read(hdfFile, '/dims/depth');

I = find(T<-10);
T(I) = NaN;

% plot
figure(1)
clf
pcolor(n(1:nSkip:end), z, T(:,1:nSkip:end))
shading interp


