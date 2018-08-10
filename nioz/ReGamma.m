clear all
close all

load apef_thorpe_wholeTM.mat

myRe = []; myGa = []; myReb = [];
for ii = 1:216*6
    ii

    I1 = find(timeSensor >= datenum(2012,10,7)+ii/(24*6) & timeSensor <= datenum(2012,10,7,0,10,0)+ii/(24*6));  

    if isempty(I1)
        continue
    end
    
    windowRho = Ritp1(:,I1);
    windowTime = timeSensor(I1);
    H = size(windowRho,1);

 
    %% Compute Re
    %Re = hst_Re('./roc12.mat',zVecReg, timeSensor1, 'zbin', 10, 'timebin', .5/24);
    Re = hst_Re('./roc12.mat',zVecReg, windowTime);
    myRe = [myRe nanmean(Re)];

    JbVec = []; 
    epsVec = [];
    for j = 1:length(windowTime)
        rhoVec_raw =  windowRho(:,j);
        [rhoVec_sort, I] = sort(rhoVec_raw);
        d = zVecReg - zVecReg(I);
        EP1 = sum(rhoVec_raw.*(max(zVecReg)-zVecReg')); 
        EP2 = sum(rhoVec_sort.*(max(zVecReg)-zVecReg'));
        rho0 = nanmean(windowRho(:,j));

        % local stratification
        N2local = g./rho0.*gradient(rhoVec_sort)./gradient(zVecReg');
        JbVec = [JbVec g./H/rho0.*(EP1-EP2).*nanmean(sqrt(N2local))];
        epsVec = [epsVec .64*rms(d).^2.*nanmean(sqrt(N2local)).^3];
        
        Reb = nanmean(epsVec)./(1e-6*nanmean(N2local));
        
    end
    myGa = [myGa nanmean(JbVec./epsVec)];
    myRe = [myRe nanmean(Re)];
    myReb = [myReb Reb];
    disp('  -> done!')

end

plot(log10(myRe), myGa, '.k')
