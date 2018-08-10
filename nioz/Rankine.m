clear
close all

skip=1;
Xvec=-100.5:skip:100.5;
%Zvec=800:skip:915;
Zvec=-100.5:skip:100.5;
R0=10;
R0=15;

[X,Z]=meshgrid(Xvec,Zvec);
R=sqrt(X.^2+Z.^2);

omega=1;
Uin=omega*R;
Uout=omega*R0^2./R;

U=Uout;U(R<R0)=Uin(R<R0);

Ux=U.*cos(atan2(X,Z));
Uz=-U.*sin(atan2(X,Z));
Ux([1:5 end-4:end],:)=0;
Ux(:,[1:5 end-4:end])=0;
Uz([1:5 end-4:end],:)=0;
Uz(:,[1:5 end-4:end])=0;

%rho=1032.5-4e-3*Z;
rho=1032.07-8e-4*Z;

rho1=rho;

nsteps=50;
dt=pi/nsteps;
for t=1:nsteps
    rho1=interp2(X,Z,rho1,X-dt*Ux,Z-dt*Uz);
end

%%
rankine.Z=Zvec;
rankine.drho=(rho1(:,ceil(end/2))-rho(:,ceil(end/2)))/((rho(end)-rho(1)));
rankine.rho = rho1;
save('rankine_pi.mat','rankine')
% save('rankine_2pi.mat','rankine')

%%
figure;
subplot(1,5,[1:4])
imagesc(Xvec,Zvec,rho1)
hold on;
contour(Xvec,Zvec,rho1,'k')
axis xy
plot([0 0],ylim,'k--')
yylim=ylim;
xlim([-90 90])
set(gca,'Xtick',[-50 0 50],'Ytick',[-50 0 50])
xlabel('x (m)');
ylabel('z (m)')
subplot(1,5,5)
plot(rho(:,ceil(end/2))-1000,Zvec);hold on;plot(rho1(:,ceil(end/2))-1000,Zvec)
ylim(yylim);
set(gca,'YtickLabel',{})
xlim([min(min(rho))-1000 max(max(rho))-1000])
xlabel('\rho-1000 (kg/m³)')
%print_figure('Rankine1','Size',[1000 1000])

figure;
imagesc(Xvec,Zvec,rho1)
hold on;
contour(Xvec,Zvec,rho1,'k')
axis xy
plot([0 0],ylim,'k--')
yylim=ylim;
xlim([-90 90])
set(gca,'Xtick',[-50 0 50],'Ytick',[-50 0 50])
xlabel('x (m)');
ylabel('z (m)')
%print_figure('Rankine2','Size',[1000 1000])


break
%%
figure;plot(rho(:,ceil(end/2)),Zvec);hold on;plot(rho1(:,ceil(end/2)),Zvec)

%%
figure;plot(rho(:,ceil(end/2))-rho1(:,ceil(end/2)),Zvec)

