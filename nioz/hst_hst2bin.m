function hst_hst2bin(inFile, outFile, dtSec)

% usage ex: 
%   iow_hst2bin('Titp_S1_p4std.mat', 'Titp_S1_p4std_1sAve.mat', 1)
%   iow_hst2bin('Titp_S1_p4std.mat', 'Titp_S1_p4std_10minAve.mat', 600)
%   iow_hst2bin('Titp_T1B_p2std.mat', 'Titp_T1B_p2std_1sAve.mat', 1)
%   iow_hst2bin('Titp_T1B_p2std.mat', 'Titp_T1B_p2std_5minAve.mat', 300)
%   hst_hst2bin('rockall_full.mat', 'Rockall_2minAve.mat', 120)
%   hst_hst2bin('Titp_T1C_p4std.mat', 'Titp_T1C_p4std_1sAve.mat', 1)
%   hst_hst2bin('Titp_T1C_p4std.mat', 'Titp_T1C_p4std_10minAve.mat', 600)
%   hst_hst2bin('Titp_T1C_p4std.mat', 'Titp_T1C_p4std_1minAve.mat', 60)


load(inFile)
dt = dtSec/86400;
timeVec = timeSensor(1)+dt/2:dt:timeSensor(end);

I = find(timeSensor>=timeVec(1) & timeSensor<=timeVec(end));
timeSensor = timeSensor(I);
Titp = Titp(:,I);

Tbin = nan(length(zVec), length(timeVec));
for i=1:length(timeVec)
    if mod(i, 1000) == 0
        disp(sprintf('%d / %d',i,length(timeVec)))
    end
    I = find(timeSensor>=timeVec(i)-dt/2 & timeSensor<=timeVec(i)+dt/2);
    Tbin(:,i) = nanmean(Titp(:,I),2);
end

save(outFile, 'Tbin', 'zVec', 'timeVec')
