function [timeVec, zVec, Titp, Traw] = hst_sensor2temp(configFile)

% function [timeVec, zVec, Titp, Traw] = hst_sensor2temp(configFile)
% 
% This function uses binary data from sensors and transfor the
% output into a matrix (t, z, T) for further analysis. Informations
% concerning the deplyment or the time window of the output is
% given into the config. file. Outputs are:
% 
% - timeVec: time vector in matlab 'mtime'
% - zVec: vertical axis (defined in the config. file)
% - Titp: Temperature matrix with missing/bad sensors interpolated
% - Traw: 
% usage ex: [timeVec, zVec, Titp, Traw] = hst_sensor2temp('ROCconfig');
%           [timeVec, zVec, Titp, Traw] = hst_sensor2temp('ROCconfig0366');


% read config file    
eval(configFile);

% load calibration files
for i = 1:size(calibFiles,1)
    fname = calibFiles(i,:);
    I = find(fname==' '); % remove empty space if there is  
                          %fname(I) = [];
    cal{i} = load(fname);
    load(fname)
    clear fname
end

% USe an existing list or generate it (see configFile)
if isempty(thermistor_list)
    thermList = hst_build_sensor_list(sensor_id, fileMinSize, '/media/Seagate1TB/NIOZ/thermistordata/ROC12/mod4c', '/media/Seagate1TB/NIOZ/thermistordata/ROC12/oldmod4');
    disp('*** A list have been generated ***') 
    disp('(look in ./sensors.list for overview)')
    disp('  -> and/or press key to continue')
    pause
else
    fid = fopen(thermistor_list);
    thermList = textscan(fid, '%s', 'delimiter', '\n'); 
end

% initialize thermistor
t=thermistor();
% This might be a OCTAVE PROB.
if size(thermList, 1) == 141
    thermList  = thermList(2:end);
end
t=append_data(t,char(thermList(1)),'nodata');

% Loop for calibration;
t_for_range=ref_date(t,REF_TIME);
range_wanted=[cdag1 cdag2];
keyboard
range=range_estimator(t_for_range,range_wanted,'yeardays');
skip=skipcal;
s=string();
for p=1:length(sensor_id)
    clear t
    t=thermistor();
    filename = char(thermList(p));
    display(filename);
    if strcmp(filename, 'missing') ~= 1 & isempty(filename) ~= 1  &  strcmp(filename, 'missing_too_short') ~= 1
        t=append_data(t,filename,'skip',skip,'range',range,'noaccel','patch');    
    end
    s=append_therm(s,t);
end
s=ref_date(s,datenum(2001,1,1));


% Calibration   (VERIF.)
% this is weak, but we suppose whether 1,2 or 3 calibfiles
if size(cal,2) == 1
    cal1 = struct2cell(cal{1});
    s=import_calib(s, cal1{1});
elseif size(cal,2) == 2
    cal1 = struct2cell(cal{1});
    cal2 = struct2cell(cal{2});
    s=import_calib(s, cal1{1}, cal2{1});
elseif size(cal,2) == 3
    cal1 = struct2cell(cal{1});
    cal2 = struct2cell(cal{2});
    cal3 = struct2cell(cal{3});
    s=import_calib(s, cal1{1}, cal2{1}, cal3{1});
else
   disp('More than 3 calib files have been found, we use 3')
   cal1 = struct2cell(cal{1});
   cal2 = struct2cell(cal{2});
   cal3 = struct2cell(cal{3});
   s=import_calib(s, cal1{1}, cal2{1}, cal3{1});
end
% replace this line:
% $$$ s=import_calib(s,S_ijkbad_2013_05_08,S_ijkbad_2013_05_16);


% Some tricks on calibration that I don't get 
a=get_value(s,'calib','start_date','id');
f=zeros(1,length(a.calib));ff=zeros(1,length(a.calib));
for k=1:length(a.calib) 
    f(k)=isempty(a.calib{k}.temperature.p);
    ff(k)=isempty(a.start_date{k}.raw.value);
end
index_nocalib=setdiff(find(f),find(ff));
a.id(index_nocalib);
s=import_calib(s,index_nocalib,cal1{1},'average'); %(VERIF)
%s=import_calib(s,index_nocalib,S_ijkbad_2013_05_08,'average');

% HERE!!!
%[s, shift]=hst_autoscale(s,'plot','treshold',thres,'method', 'polyval','degree',degr, 'shift'); %Test
                                                                                                

% auto-scale
%s=auto_scale(s,'plot','treshold',thres,'method','polyval','degree',degr);%bathcalib
[s, shift]=auto_scale(s,'plot','treshold',thres,'method', 'polyval','degree',degr, 'shift');%bathcalib

% Autoscale 2 to correct few sensors
I = find(abs(shift)>.02); % this could be optimized
%s2=auto_scale2(s,I,'n_iter',1);                          
s2=auto_scale2(s,[13 15 65 85 110],'n_iter',1);
s3=auto_scale2(s2,[13 15],'n_iter',1); % test, not working
                                       % really... - FC


% Here you can by-pass auto_scale2:  % -FC
s_calib=remove_data(s2);
%s_calib=remove_data(s);


% Loop to build the string (shorter period):
s=string();
range_wanted=[dag1 dag2]; %read short period
range=range_estimator(t_for_range,range_wanted,'yeardays')
skip=skipdata;
s=string();
for p=1:length(sensor_id)
    clear t
    t=thermistor();
    filename = char(thermList(p));
    display(filename);
    if strcmp(filename, 'missing') ~= 1 & isempty(filename) ~= 1  &  strcmp(filename, 'missing_too_short') ~= 1
        t=append_data(t,filename,'skip',skip,'range',range,'noaccel','patch');    
    end
    s=append_therm(s,t);
end
s=ref_date(s,datenum(2001,1,1));
s=import_calib(s,s_calib);

% Build the matrix
[tijd,Torig]=string2mat(s,'calib','yearday');
TT = Torig; %[tijd,TT]=string2mat(s,'calib','yearday');
[mT,nT]=size(TT);
timeVec = tijd + datenum(YYYY, 1, 1, 0, 0, 0);

% Test coment
% $$$ sa=35.3*ones(mT,nT);
% $$$ dplct10=gsw_CT_from_t(sa,TT,dz*1.01-10.1325);
% $$$ TT=dplct10;
% $$$ nf1=1;
% $$$ nf2=nT;

% midT=mean(TT(:,2:nf2-1)');
% adiab=1.86876e-4;
% adiabp=adiab*[1:103]+midT(1);
% stabtrend=(midT(mT)-adiabp(mT)-(midT(1)-adiabp(1)))/(mT-1);
% figure(1);clf;plot(midT,dz,adiabp,dz,'r')
% pause

% Interpolation between missing sensors (weak, but works here)
tempVec = nanmean(Torig,2);
if isempty(Iremove)
    % smooth fit
    I = find(~isnan(tempVec));
    p = polyfit(I, tempVec(I), 3);
    x = [1:length(tempVec)]';
    tempFit = p(1)*x.^3 + p(2)*x.^2 + p(3)*x + p(4);
    
    % remove outliers
    sd = nanstd(abs(tempVec - tempFit));
    I = find(abs(tempVec-tempFit) > 3*sd);
    tempVec(I) = NaN;

    % Show missing values
    I = find(isnan(tempVec));
    disp([' !! No pts to remove was suggested. We suggest removing ' ...
          'these sensors: !!'])
    sprintf( '%d ', I)
else
    tempVec(Iremove) = NaN; % use existing remove list
end    

% the intperolation
I = find(~isnan(tempVec)==1);
for i = 1:size(TT, 2)
    TT(:,i) = interp1(I, TT(I,i), 1:size(TT,1));
end

if interactive
    figure(1)
    clf
    imagesc(TT)
    bool = 1;
    while bool
        R = input('Does Figure 1 looks okay? and do you want to continue (y/n)', 's');
        I = find(R==' '); % remove white space
        R(I) = [];
        if strcmp(R, 'n') == 1
            disp( '-> ABORT NOW!')
            return         
        elseif strcmp(R, 'y') == 1
            bool = 0;
        else
            disp(' please answer y or n')
        end
    end
end
   
Titp = TT;
Traw = TT;
Traw(Iremove, :) = NaN;

figure(1)
clf
contourf(timeVec, zVec, Titp, 100, 'linestyle', 'none')
set(gca,'ydir', 'reverse')
datetick
xlabel(datestr(timeVec(1), 1))
ylabel('Depth (m)')
xlim([min(timeVec) max(timeVec)])

figname = sprintf('./figures/temp_%s.png', datestr(timeVec(1), 29));
print(gcf, '-dpng', '-r100', figname)
