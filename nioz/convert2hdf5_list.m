function convert2hdf5_list(filelist, ext)
    
    
% load *.P files names (file in which are recorded .P files)
fid = fopen(filelist);
C = textscan(fid, '%s', 'delimiter', '\n');
files = char(C{1});


noFiles = size(files, 1);

for i = 1:noFiles

    fname = files(i, :);
    I = find(fname==' ');   
    fname(I) = [];

    convert2hdf5_fc(fname,ext)
    
end
