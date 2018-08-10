function spm_moorings(listfile, tlim, zmin)

% usage ex:
% spm_moorings('MIQ60.list', [datenum(2017,9,6,18,0,0), datenum(2017,9,13,17,0,0)], -65)
% spm_moorings('MIQ40.list', [datenum(2017,8,30,11,0,0), datenum(2017,9,13,17,0,0)], -38)
% spm_moorings('MIQ40.list', [datenum(2017,9,5,0,0,0), datenum(2017,9,8,0,0,0)],-38)

% Time ,limits
t0 = tlim(1);
tf = tlim(2);
dt=1/(24*60); % 
t=[t0:dt:tf]';

% open list
fid = fopen(listfile);
C = textscan(fid, '%s', 'delimiter', '\n');
files = char(C{1});
no_files = size(files, 1); %number of *.P files 

% loop on files
Mast2D = struct();
for i=1:no_files
    fname = files(i, :);
    [jj,mm,yyyy,HH,MM,SS,bid,temp0,pression,tension,bid1,temp1bid] =textread(fname,'Date: %f/%f/%f %f:%f:%f %s %f P: %f Tension: %f %s %f' ,'headerlines',15); 
    Mast2D(i).date=datenum(yyyy,mm,jj,HH,MM,SS);
    Mast2D(i).temp=temp0;
    Mast2D(i).pres=pression;
end
% pressure now in kPa

%% ---- Data check plots ---- %% 
% Temperature timeseries
% $$$ figure(1);clf;hold on
% $$$ set(1,'position',[10 10 1500 500])
% $$$ for nsta=1:no_files
% $$$ 
% $$$     plot(Mast2D(nsta).date(:),Mast2D(nsta).temp,'k','linewidth',2)
% $$$ 
% $$$     tlabel
% $$$     title(num2str(nsta))
% $$$     grid on
% $$$ 
% $$$ end

% Pressure timeseries
figure(2);clf;hold on
set(2,'position',[10 10 1500 500])

for nsta=1:no_files

    if mean(Mast2D(nsta).pres) > 10000
        Mast2D(nsta).pres=NaN(size(Mast2D(nsta).temp));
    end
    offset=mean(Mast2D(nsta).pres(1:10));
    Mast2D(nsta).pres=-(Mast2D(nsta).pres-offset)*10; %<-- conversion to dbar
    plot(Mast2D(nsta).date(:),Mast2D(nsta).pres,'k','linewidth',2)

    tlabel
    title(num2str(nsta))
    grid on

end

% $$$ % Color scatter plot
% $$$ figure(4);clf;hold on
% $$$ set(4,'position',[10 10 1500 500])
% $$$ for n=1:no_files
% $$$     scatter(Mast2D(n).date(:),Mast2D(n).pres(:),30,Mast2D(n).temp(:),'o','filled')
% $$$ end
% $$$ colorbar
% $$$ tlabel
% $$$ caxis([2 16])
% $$$ %xlim([ tdeb tfin])
% $$$ grid on


% pour commencer, on vire les capteurs de T° seuls
for n=1:no_files
    Mast2D(n).tempi=interp1(Mast2D(n).date,Mast2D(n).temp,t);
    if sum(isnan(Mast2D(n).pres))== 0  
        Mast2D(n).presi=interp1(Mast2D(n).date,Mast2D(n).pres,t);
    else
        Mast2D(n).presi=NaN(size(t));
    end
end


% on verifie avec un scatter
% $$$ figure(5);clf;hold on
% $$$ set(5,'position',[10 10 1500 500])
% $$$ for n=1:no_files
% $$$     scatter(t,Mast2D(n).presi(:),30,Mast2D(n).tempi(:),'o','filled')
% $$$ end
% $$$ colorbar
% $$$ tlabel
% $$$ caxis([ 2 17])
% $$$ %xlim([ tdeb tfin])
% $$$ grid on


% $$$ % on verifie que l'interpolation temporelle s'est bien passée
% $$$ figure;hold on
% $$$ for n=1:no_files
% $$$     plot(t,Mast2D(n).tempi)
% $$$ end
% $$$ grid on
% $$$ title('apres interpolation temporelle')

nbrec=size(t,1);
mmax=nbrec
% on interpole sur la verticale
z=[zmin:2:0]';
nbniv=size(z,1);

temp_inter=zeros(mmax,size(z,1));
m=0;
for n=1:nbrec
    m=m+1;
    for k=1:no_files
        zz(k)=Mast2D(k).presi(n);
        temp1(k)=Mast2D(k).tempi(n);
    end
    date_inter(m)=t(n);
    temp_inter(m,:)=interp1(zz,temp1,z');
end


figure(100);clf;hold on
set(100,'position',[10 10 1500 500])

%imagesc(date_inter,z',temp_inter')
contourf(date_inter,z',temp_inter',40,'linestyle','none')
[C, H] = contour(date_inter,z',temp_inter',[13,13]);
plot(C(1,:), C(2,:), 'k')

%tlabel('x',' dd/mm HH')
tlabel
grid on
set(gca,'layer','top','Linewidth',1.1,'gridlinestyle','--','fontsize',12)
h=colorbar('location','southoutside');

set(0,'DefaultFigurePaperPositionMode','auto')
%  title([ fname ' Cross shore raw baroclinic velocity
%  '],'fontsize',14,'interpreter','none')
axis([t0 tf zmin 0])
caxis([2 16])
ylabel('Profondeur(m)','fontsize',16)
title(' Rade Miquelon','fontsize',18)
%axis( min(date_inter(1)) max(date_inter(end)) 0 30])
repfig='./';
set(0,'DefaultFigurePaperPositionMode','auto')
ficfig=fullfile(repfig,' rade_miquelon_60m');
print(100,'-dpng',ficfig)  


C(:, 1:2) = [];
I = find(C(2, :)>40);
C(:, I) = [];


keyboard
figure(13)
clf
%plot(C(1,:), diff(C(2,:)), 'k')
plot(C(1,:), C(2,:)-runmean(C(2,:), 60), 'k')
%ylim([-3, 3]) 
xlim([datenum(201)])
