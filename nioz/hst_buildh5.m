function hst_buildh5(configFile)

% usage ex: hst_buildh5('ROCconfig')

% VERIF means that the portion of script has been verified that it
% was not inducing error...
    
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
    C = textscan(fid, '%s', 'delimiter', '\n'); 
end

% initialize thermistor
t=thermistor();
t=append_data(t,char(thermList(1)),'nodata');

% Loop for calibration;
t_for_range=ref_date(t,REF_TIME);
range_wanted=[cdag1 cdag2];
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




keyboard




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

