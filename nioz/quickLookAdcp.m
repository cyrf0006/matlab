% Assume ADCP structure exists and is loaded in memory
% Ex: for Albex lander 2013

zVec = ADCP.config.ranges;
timeVec = ADCP.mtime;
E = ADCP.east_vel;
N = ADCP.north_vel;
W = ADCP.vert_vel;
Int1 = squeeze(ADCP.intens(:,1,:));
Int2 = squeeze(ADCP.intens(:,2,:));
Int3 = squeeze(ADCP.intens(:,3,:));
Int4 = squeeze(ADCP.intens(:,4,:));


%% Correct for magnetic declination!
% (declination is 8.415 deg. West)
[E, N] = rotate_vecd(E, N, 8.4150);

%% reduce time 
t0 = datenum(2013, 10, 8, 19, 30, 0);
tf = datenum(2013, 10, 9, 16, 30, 0);

I = find(timeVec>=t0 & timeVec<=tf);

%imagesc(timeVec(I), zVec, Int(:,I));
%imagesc(timeVec(I), zVec, Int(:,I));

IntMean1 = nanmean(Int1(:,I),2);
IntMean2 = nanmean(Int2(:,I),2);
IntMean3 = nanmean(Int3(:,I),2);
IntMean4 = nanmean(Int4(:,I),2);
for i = 1:size(Int1,1)
    Int1(i,:) = Int1(i,:)-IntMean1(i);
    Int2(i,:) = Int2(i,:)-IntMean2(i);
    Int3(i,:) = Int3(i,:)-IntMean3(i);
    Int4(i,:) = Int4(i,:)-IntMean4(i);
end

[U,V] = rotate_vecd(E,N,-30);
load BR_symetric
figure(1)
clf
subplot(411)
plot(timeVec(I), ADCP.temperature(I))
datetick
xlim([t0 tf])
set(gca, 'xticklabel', [])
ylabel('T(^{\circ}C)')
POS1 = get(gca, 'pos');

subplot(4,1,[2 4])
imagesc(timeVec(I), zVec, V(:,I));
set(gca, 'ydir', 'normal')
datetick
xlim([t0 tf])
ylim([0 25])
colorbar
colormap(BR_symetric)
POS2 = get(gca, 'pos');
POS2(3) = POS1(3);
set(gca, 'pos', POS2)

figure(2)
clf
subplot(411)
plot(timeVec(I), ADCP.temperature(I))
datetick
xlim([t0 tf])
set(gca, 'xticklabel', [])
ylabel('T(^{\circ}C)')
POS1 = get(gca, 'pos');

subplot(4,1,[2 4])
imagesc(timeVec(I), zVec, Int3(:,I))
set(gca, 'ydir', 'normal')
%caxis([0 50])
datetick
xlim([t0 tf])
ylim([0 25])
colorbar
colormap(BR_symetric)
POS2 = get(gca, 'pos');
POS2(3) = POS1(3);
set(gca, 'pos', POS2)

figure(3)
clf
plot(nanmean(V(:,I(1:320)), 2), zVec)
hold on
plot(nanmean(V(:,I(320:end)), 2), zVec, 'r')
ylim([0 25])


figure(4)
clf
subplot(411)
plot(timeVec(I), ADCP.temperature(I))
datetick
xlim([t0 tf])
set(gca, 'xticklabel', [])
ylabel('T(^{\circ}C)')
POS1 = get(gca, 'pos');

subplot(4,1,[2 4])
imagesc(timeVec(I), zVec, W(:,I));
set(gca, 'ydir', 'normal')
datetick
xlim([t0 tf])
ylim([0 25])
caxis([-.05 .05])
colorbar
colormap(BR_symetric)
POS2 = get(gca, 'pos');
POS2(3) = POS1(3);
set(gca, 'pos', POS2)


%% Shear
figure(5)
clf

I = find(timeVec>=t0 & timeVec<=tf);
J = find(zVec<=25);
U = U(J,I);
V = V(J,I);
zVec = zVec(J);
timeVec = timeVec(I);
tempVec = ADCP.temperature(I);
% itp vel
for i = 1:length(zVec);
    uVec = U(i,:);
    I = find(~isnan(uVec));
    U(i,:) = interp1(timeVec(I), uVec(I), timeVec);
    
    vVec = V(i,:);
    I = find(~isnan(vVec));
    V(i,:) = interp1(timeVec(I), vVec(I), timeVec);
end

[FX, duz] = gradient(U);
[FX, dvz] = gradient(V);



dudz = nan(size(duz));
dvdz = nan(size(duz));
for i = 1:length(timeVec)
    dudz(:,i) = duz(:,i)./gradient(zVec,1,1);
    dvdz(:,i) = dvz(:,i)./gradient(zVec,1,1);
end
S2 = dudz.^2 + dvdz.^2;


subplot(411)
plot(timeVec, tempVec)
datetick
xlim([t0 tf])
set(gca, 'xticklabel', [])
ylabel('T(^{\circ}C)')
POS1 = get(gca, 'pos');

subplot(4,1,[2 4])
imagesc(timeVec, zVec, log10(S2));
set(gca, 'ydir', 'normal')
datetick
xlim([t0 tf])
ylim([0 25])
%caxis([-.05 .05])
colorbar
colormap(BR_symetric)
POS2 = get(gca, 'pos');
POS2(3) = POS1(3);
set(gca, 'pos', POS2)

