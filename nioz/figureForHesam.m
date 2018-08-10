clear all
close all

load Hesam.mat

dt = 10/60/24;
timeVec = round(timeSensor(1)*24)/24:dt/2:timeSensor(end);
Ri = hst_Ri('./roc12.mat',N2bg,zVecReg, timeVec);
RiVecRaw = nanmean(Ri.Ri, 1);
Re = hst_Re('./roc12.mat',zVecReg, maxdVec, timeVec);

TVec = nanmean(T,1);
nu = visc35(TVec);
Reb = epsVec./(nu.*nanmean(N2bg,1));

GammaVec = nan(size(timeVec));
ReVec = nan(size(timeVec));
RiVec = nan(size(timeVec));
for i = 1:length(timeVec)   
    I = find(timeSensor>=timeVec(i)-dt/2 & timeSensor<timeVec(i)+dt/2);
    GammaVec(i) = nanmean(JbVec(I))./nanmean(epsVec(I));
    RebVec(i) = nanmean(Reb(I));

    I = find(Ri.timeVec>=timeVec(i)-dt/2 & Ri.timeVec<timeVec(i)+dt/2); 
    ReVec(i) = nanmean(Re(I));
    RiVec(i) = nanmean(RiVecRaw(I));        
end

% $$$ size(GammaVec)
% $$$ size(ReVec)
% $$$ size(RiVec)

figure(13)
clf
subplot(131)
semilogx(RiVec, GammaVec, '.k')
xlabel('Ri')
ylabel('\gamma')
xlim([1e-2 5e5])
set(gca, 'xgrid', 'on', 'ygrid', 'on')

subplot(132)
semilogx(ReVec, GammaVec, '.k')
xlabel('Re')
ylabel('\gamma')
set(gca, 'xgrid', 'on', 'ygrid', 'on')

subplot(133)
semilogx(RebVec, GammaVec, '.k')
xlabel('Re_b')
ylabel('\gamma')
xlim([5e2 8e5])
set(gca, 'xgrid', 'on', 'ygrid', 'on')
