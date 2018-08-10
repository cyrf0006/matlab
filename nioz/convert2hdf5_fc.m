function convert2hdf5_fc(filename_read,type,varargin)
% convert2hdf5(filename_read,type) will rewrite the binary data from the
% .dec file filename_read in a hdf5 format.
% type ='v3','v4','v4b', 'v4c' or 'v5' is used to precise the
% type of file which is on input.
% if filename_read is a directory, the function will convert all the .dec
% files found in the directory.
% convert2hdf5(filename_read,type,'force') will force the overwriting
% of the files.
% convert2hdf5(filename_read,type,'ignore_nsamples') will ignore the
% information on the number of data points from the header but compute the
% number of samples from the size of the file.
%
% v4 : very first v4 output
% v4b : sensors with compass for km3net
% v4c : current v4 output
% v4d : with 8 bytes for prog version before id
% v5b : new version, for spain for instance 
% ts4 : 


%% browse options
p=1;
force=0;
ignore_nsamples=0;
while p<=length(varargin)
    arg=varargin{p};
    if isstr(arg)
        switch arg
            case {'force'}
                force=1;
                p=p+1;
            case {'ignore_nsamples'}
                ignore_nsamples=1;
                p=p+1;
            otherwise
                fprintf('argument %s not understood\n',arg)
                return
        end
    end
end


%% check the file type is supported
if ~ismember(type,{'v3','v3b','v3c','v4','v4b','v4c','v4d','v5','v5b','ts4'})
    display(['Thermistor type ' type ' unknown; use v3, v3b, v3c v4, v4b, v4c, v4d, v5, v5b or ts4 instead !'])
    return
end

if ~isdir(filename_read)
%     switch type
%         case 'ts4'
            filename_write=regexprep(filename_read,'(.*\.)\w+$',['$1' 'hdf']);
%         otherwise
%             filename_write=[filename_read '.' type];
%     end
    if exist(filename_write)
        if force
            delete(filename_write)
        else
            display(['file ' filename_write ' already exists'])
            return
        end 
    end
end

if ~exist(filename_read)
    display(['file ' filename_read ' does not exists'])
    return
end

%% in case the input file is a folder, browse recursively the folder
if isdir(filename_read)
    f=dir(filename_read);
    for i=1:length(f)
        switch f(i).name
            case {'.','..'}
                continue
        end
        if isdir(f(i).name)            
            convert2hdf5(fullfile(filename_read,f(i).name),type,varargin{:});
            return
        else
            [p,n,e] = fileparts(f(i).name);
            switch lower(e)
                case {'.dec','.ts4'}
                    convert2hdf5(fullfile(filename_read,f(i).name),type,varargin{:});
            end
        end
    end
    return
end

switch type
    case 'ts4'
        tic;
        display(['Converting file ' filename_read ' with ts4dump.py']);            
        system([fullfile(fileparts(which('convert2hdf5')),'ts4dump.py')...
            ' -chdf -f ' filename_read]);
        display(['File ' filename_read ' converted in ' num2str(toc) 's.']);            
        return
end

%% open the input and output files
fid_read=fopen(filename_read);
% fid_write=fopen(filename_write,'w+');
keyboard
%h5create(filename_write,'/raw_data/temperature',Inf,'ChunkSize',8192);
hdf5write(filename_write,'/raw_data/temperature',Inf,'ChunkSize',8192);


%% start conversion
tic;

%% search for the numb of samples in the footer and the eof position
switch type
    case {'v3','v3b','v3c','v4','v4b','v4c','v4d'}
        fseek(fid_read,-8,'eof');
        n_samples=fread(fid_read,1,'uint32');
    case 'v5'
        fseek(fid_read,-8,'eof');
        n_samples=fread(fid_read,1,'uint32');
        n_samples=round(n_samples/2);
    case 'v5b'
        fseek(fid_read,-4,'eof');
        n_samples=fread(fid_read,1,'uint32');
%         n_samples=round(n_samples/2);
end

pos_end=ftell(fid_read);




%% read data from header
fseek(fid_read,0,'bof');
% fwrite(fid_write,hex2dec('FFFF'),'uint16');%header tag 'FFFF'
switch type
    case 'v3'
        fwrite(fid_write,fread(fid_read,30));
    case 'v3b'
        fwrite(fid_write,fread(fid_read,38));
    case 'v3c'
        fseek(fid_read,8,'cof');     
        fwrite(fid_write,fread(fid_read,30));
    case {'v4','v4b'}
        fwrite(fid_write,fread(fid_read,44));
        fseek(fid_read,42,'bof');
        n_skip_accel=fread(fid_read,1,'uint16');
        accel_count=n_skip_accel-1;
    case {'v4c'}
        fwrite(fid_write,fread(fid_read,36));
        fseek(fid_read,12,'bof');
        n_skip_accel=fread(fid_read,1,'uint16');
        accel_count=n_skip_accel-1;
    case {'v4d'}
        h5writeatt(filename_write,'/','FW_VersionDate_s',fread(fid_read,1,'uint32'));
        h5writeatt(filename_write,'/','FW_VersionMajor',fread(fid_read,1,'uint8'));
        h5writeatt(filename_write,'/','FW_VersionMinor',fread(fid_read,1,'uint8'));
        h5writeatt(filename_write,'/','UnitType',fread(fid_read,1,'uint16'));
        h5writeatt(filename_write,'/','MyID',fread(fid_read,1,'uint16'));
        h5writeatt(filename_write,'/','AliveTO',fread(fid_read,1,'uint8'));
        h5writeatt(filename_write,'/','ComTestTime_s',fread(fid_read,1,'uint8'));
        h5writeatt(filename_write,'/','IMsignalTO_s',fread(fid_read,1,'uint8'));
        h5writeatt(filename_write,'/','ConfirmT_s',fread(fid_read,1,'uint8'));
        h5writeatt(filename_write,'/','HF_StartTime_ms',fread(fid_read,1,'uint8'));
        h5writeatt(filename_write,'/','LF_StartTime_ms',fread(fid_read,1,'uint8'));
        h5writeatt(filename_write,'/','SD_PwrCycling',fread(fid_read,1,'uint8'));
        h5writeatt(filename_write,'/','ADC_SamplesPMeas_2comp',fread(fid_read,1,'uint8'));
        h5writeatt(filename_write,'/','ACCEL_StartTime_ms',fread(fid_read,1,'uint8'));
        h5writeatt(filename_write,'/','Dummy_0',fread(fid_read,1,'uint8'));
        h5writeatt(filename_write,'/','Accel_SkipCnt',fread(fid_read,1,'uint16'));
        h5writeatt(filename_write,'/','BAT_MinVoltage_mv',fread(fid_read,1,'uint16'));
        h5writeatt(filename_write,'/','SampleInt_ms',fread(fid_read,1,'uint16'));
        h5writeatt(filename_write,'/','SampleTime_ms',fread(fid_read,1,'uint16'));
        h5writeatt(filename_write,'/','ShortSyncInt_s',fread(fid_read,1,'uint16'));
        h5writeatt(filename_write,'/','ShortSyncCnt',fread(fid_read,1,'uint16'));
        h5writeatt(filename_write,'/','LongSyncInt_s',fread(fid_read,1,'uint16'));
        h5writeatt(filename_write,'/','CruiseID',fread(fid_read,1,'uint16'));
        h5writeatt(filename_write,'/','Station',fread(fid_read,1,'uint16'));
        h5writeatt(filename_write,'/','NewTime_s',fread(fid_read,1,'uint32'));
        h5writeatt(filename_write,'/','ResendTO',fread(fid_read,1,'uint16'));

        n_skip_accel=h5readatt(filename_write,'/','Accel_SkipCnt');
        accel_count=n_skip_accel-1;

    case {'v5'}
        fwrite(fid_write,fread(fid_read,42));
        fseek(fid_read,16,'bof');
        n_skip_accel=fread(fid_read,1,'uint16');
        accel_count=n_skip_accel-1;
    case {'v5b'}
        fseek(fid_read,8,'cof')
        fwrite(fid_write,fread(fid_read,42));
        fseek(fid_read,8+16,'bof');
        n_skip_accel=fread(fid_read,1,'uint16');
        accel_count=n_skip_accel-1;
end


switch type
    case 'v3'
        n_samples_guess=floor((pos_end-108-30)/4);
    case 'v3b'
        n_samples_guess=floor((pos_end-108-38)/4);
    case 'v3c'
        n_samples_guess=floor((pos_end-108-30-8)/4);
    case {'v4'}
        if n_skip_accel>0
            n_samples_guess=floor((pos_end-108-44)/10);
        else
            n_samples_guess=floor((pos_end-108-44)/4);
        end
    case {'v4b'}
        if n_skip_accel>0
        n_samples_guess=floor((pos_end-108-44)/16);
        else
        n_samples_guess=floor((pos_end-108-44)/4);
        end
    case {'v4c'}
        if n_skip_accel>0
        n_samples_guess=floor((pos_end-108-36)/10);
        else
        n_samples_guess=floor((pos_end-108-36)/4);
        end
    case {'v4d'}
        if n_skip_accel>0
            n_samples_guess=floor((pos_end-108-36-8)/10);
        else
            n_samples_guess=floor((pos_end-108-36-8)/4);
        end
    case {'v5'}
        if n_skip_accel>0
            n_samples_guess=floor((pos_end-108-42)/14);
        else
            n_samples_guess=floor((pos_end-108-42)/8);
        end
    case {'v5b'}
        if n_skip_accel>0
            n_samples_guess=floor((pos_end-104-50)/14);
        else
            n_samples_guess=floor((pos_end-104-50)/8);
        end
end

if n_samples_guess<0
    n_samples_guess=0;
end

fprintf('File size (bytes) : %d\n',pos_end)
fprintf('Number of samples from footer : %d\n',n_samples)
fprintf('Number of samples from filesize : %d\n',n_samples_guess)

if ignore_nsamples
    n_samples=n_samples_guess;
    fprintf('  ignore_nsamples=1 : value used for n_samples = %d\n\n',n_samples)
elseif abs(n_samples-n_samples_guess)>100
    n_samples=n_samples_guess;
    fprintf('  Number of samples mismatch : value used for n_samples = %d\n\n',n_samples)    
end


pos_accel=[];
pos_cond=[];
pos_msg=[];



%% prepare the file for writing and the position vectors (fill with 0)
switch type
%     case 'v3'
%         fseek(fid_read,30,'bof');
%         fwrite(fid_write,zeros(n_samples,1),'uint32');% fill with 0
%         fseek(fid_write,-4*n_samples,'cof');% prepare for writing
%     case 'v3b'
%         fseek(fid_read,38,'bof');
%         fwrite(fid_write,zeros(n_samples,1),'uint32');% fill with 0
%         fseek(fid_write,-4*n_samples,'cof');% prepare for writing
%     case 'v3c'
%         fseek(fid_read,30+8,'bof');
%         fwrite(fid_write,zeros(n_samples,1),'uint32');% fill with 0
%         fseek(fid_write,-4*n_samples,'cof');% prepare for writing
%     case {'v4'}
%         fseek(fid_read,44,'bof');
%         dummy=zeros(n_samples,1,'uint32');
%         fwrite(fid_write,dummy,'uint32');% fill with 0
%         clear dummy;
% %         n_bytes=n_samples*4+floor(n_samples/max(n_skip_accel,1))*6;
%         if n_skip_accel>0
%             n_accel=0;
%             dummy=zeros(3*ceil(n_samples/n_skip_accel)+1,1,'uint16');
%             fwrite(fid_write,dummy,'uint16');% fill with 0
%             clear dummy ;
%             fseek(fid_write,-2*(3*ceil(n_samples/n_skip_accel)+1),'cof');% prepare for writing
%             pos_accel=zeros(ceil(n_samples/n_skip_accel)+1,1,'uint32');
%         end
%         fseek(fid_write,-4*n_samples,'cof');% prepare for writing
%     case {'v4b'}
%         fseek(fid_read,44,'bof');
% %         n_bytes=n_samples*4+floor(n_samples/max(n_skip_accel,1))*12;
%         fwrite(fid_write,zeros(n_samples,1),'uint32');% fill with 0
%         if n_skip_accel>0
%             n_accel=0;
%             pos_accel=zeros(ceil(n_samples/n_skip_accel)+1,1,'uint16');
%             fwrite(fid_write,zeros(6*ceil(n_samples/n_skip_accel)+1,1),'uint16');% fill with 0
%             fseek(fid_write,-2*(6*ceil(n_samples/n_skip_accel)+1),'cof');% prepare for writing
%         end
%         fseek(fid_write,-4*n_samples,'cof');% prepare for writing
%     case {'v4c'}
%         fseek(fid_read,36,'bof');
% %         n_bytes=n_samples*4+floor(n_samples/max(n_skip_accel,1))*6;
%         fwrite(fid_write,zeros(n_samples,1),'uint32');% fill with 0
%         if n_skip_accel>0
%             n_accel=0;
%             pos_accel=zeros(ceil(n_samples/n_skip_accel)+1,1,'uint32');
%             fwrite(fid_write,zeros(3*ceil(n_samples/n_skip_accel)+1,1),'uint16');% fill with 0
%             fseek(fid_write,-2*(3*ceil(n_samples/n_skip_accel)+1),'cof');% prepare for writing
%         end
%         fseek(fid_write,-4*n_samples,'cof');% prepare for writing
     case {'v4d'}
         fseek(fid_read,36+8,'bof');
% %         n_bytes=n_samples*4+floor(n_samples/max(n_skip_accel,1))*6;
%         fwrite(fid_write,zeros(n_samples,1),'uint32');% fill with 0
%         if n_skip_accel>0
%             n_accel=0;
%             pos_accel=zeros(ceil(n_samples/n_skip_accel)+1,1,'uint32');
%             fwrite(fid_write,zeros(3*ceil(n_samples/n_skip_accel)+1,1),'uint16');% fill with 0
%             fseek(fid_write,-2*(3*ceil(n_samples/n_skip_accel)+1),'cof');% prepare for writing
%         end
%         fseek(fid_write,-4*n_samples,'cof');% prepare for writing
%     case 'v5'
%         fseek(fid_read,42,'bof');
%         fwrite(fid_write,zeros(n_samples,1),'uint32');% fill with 0
%         fwrite(fid_write,zeros(n_samples,1),'uint32');% fill with 0
%         n_cond=0;
%         pos_cond=zeros(n_samples,1,'uint32');
%         if n_skip_accel>0
%             n_accel=0;
%             pos_accel=zeros(ceil(n_samples/max(n_skip_accel,1))+1,1,'uint32');
%             fwrite(fid_write,zeros(3*(ceil(n_samples/n_skip_accel)+1),1),'uint16');% fill with 0
%             fseek(fid_write,-2*3*(ceil(n_samples/n_skip_accel)+1),'cof');% prepare for writing
%         end
%         fseek(fid_write,-4*n_samples-4*n_samples,'cof');% prepare for writing
%     case 'v5b'
%         fseek(fid_read,50,'bof');
%         fwrite(fid_write,zeros(n_samples,1),'uint32');% fill with 0
%         fwrite(fid_write,zeros(n_samples,1),'uint32');% fill with 0
%         n_cond=0;
%         pos_cond=zeros(n_samples,1,'uint32');
%         if n_skip_accel>0
%             n_accel=0;
%             pos_accel=zeros(ceil(n_samples/max(n_skip_accel,1))+1,1,'uint32');
%             fwrite(fid_write,zeros(ceil(n_samples/n_skip_accel)+1,1),'uint16');% fill with 0
%             fwrite(fid_write,zeros(ceil(n_samples/n_skip_accel)+1,1),'uint16');% fill with 0
%             fwrite(fid_write,zeros(ceil(n_samples/n_skip_accel)+1,1),'uint16');% fill with 0
%             fseek(fid_write,-2*3*(ceil(n_samples/n_skip_accel)+1),'cof');% prepare for writing
%         end
%         fseek(fid_write,-4*n_samples-4*n_samples,'cof');% prepare for writing
end
pos_msg=[];

%% write T data and find error message positions
p=0;
t1=tic;

% fwrite(fid_write,hex2dec('FFF0'),'uint16');% t_data tag 'FFF0'
n_samples
dt1=zeros(n_samples,1);

while ftell(fid_read)<pos_end-109
    dt1(p+1)=fread(fid_read,1,'uint32');
    
%     fwrite(fid_write,dt1,'uint32');
%     h5write(filename_write,'/raw_data/temperature',bitand(dt1,hex2dec('fffffff')),p+1,1);
    
    switch type
        case {'v5','v5b'}
            n_cond=n_cond+1;
            pos_cond(n_cond)=ftell(fid_read);
            fseek(fid_read,4,'cof');
    end
    if bitand(dt1(p+1),hex2dec('80000000'))
        pos_msg(length(pos_msg)+1)=ftell(fid_read);
        fseek(fid_read,4,'cof');
    end

    
    switch type
        case {'v4','v4b','v5','v4c','v4d','v5b'}
            if accel_count>0
                accel_count=accel_count-1;
            else if accel_count==0
                    accel_count=n_skip_accel-1;
                    n_accel=n_accel+1;
                    pos_accel(n_accel)=ftell(fid_read);
                    switch type
                        case {'v4','v5','v4c','v4d','v5b'}
                            fseek(fid_read,6,'cof');
                        case {'v4b'}
                            fseek(fid_read,12,'cof');
                    end
                end
            end
    end
    if p>n_samples
        display(p)
        ftell(fid_read)
    end
    if mod(p,round(n_samples/100))==1
        display(['Temperature for file ' filename_read ' parsed at ' num2str(p/n_samples*100,'%02.2f') '% in ' num2str(toc) 's.']);        
    end
    p=p+1;
end



tic;
h5write(filename_write,'/raw_data/temperature',bitand(dt1,hex2dec('fffffff')),1,length(dt1));
display(['Temperature written to hdf5 file ' filename_write ' in ' num2str(toc) 's.']);  

switch type
    case {'v4','v4b','v4c','v4d'}
        if n_skip_accel>0
            pos_accel(n_accel+1:end)=[];
        end
    case {'v5','v5b'}
        if n_skip_accel>0
            pos_accel(n_accel+1:end)=[];
        end
        pos_cond(n_cond+1:end)=[];
end

time_temp_read=toc(t1);
tic;


%% write acceleration data
switch type
    case {}%'v4','v5','v4b','v4c','v4d','v5b'}
        if n_skip_accel>0
        fwrite(fid_write,hex2dec('FF0F'),'uint16');%accel_data tag 'FF0F'
        for k=1:length(pos_accel)
            fseek(fid_read,double(pos_accel(k)),'bof');
            fwrite(fid_write,fread(fid_read,3,'uint16'),'uint16');
            if mod(k,round(n_samples/n_skip_accel/100))==1
                display(['Acceleration for file ' filename_read ' parsed at ' num2str(k/n_skip_accel/n_samples*100,'%02.2f') '% in ' num2str(toc) 's.']);
            end
        end
        end
%        fwrite(fid_write,hex2dec('FF0F'),'uint16');%accel_data tag 'FF0F'
%         accel_stamp=zeros(6000,1,'uint16');
%         cursor=1;
%         for k=1:length(pos_accel)
%             fseek(fid_read,double(pos_accel(k)),'bof');
%             accel_stamp(cursor:cursor+2)=fread(fid_read,3,'uint16');
%             cursor=cursor+3;
%             if cursor>6000
%                 fwrite(fid_write,accel_stamp,'uint16');
%                 cursor=1;
%             end
%             if mod(k,round(n_samples/n_skip_accel/100))==1
%                 display(['Acceleration for file ' filename_read ' parsed at ' num2str(k/n_skip_accel/n_samples*100,'%02.2f') '% in ' num2str(toc) 's.']);
%             end
%             if cursor>1
%                 fwrite(fid_write,accel_stamp(1:cursor-1),'uint16');
%             end
%         end
        
end

time_accel_read=toc;
tic;

%% write compass data
switch type
    case {}%'v4b'}
        fwrite(fid_write,hex2dec('Fts4'),'uint16');%compass_data tag 'Fts4'
        for k=1:length(pos_accel)
            fseek(fid_read,double(pos_accel(k))+6,'bof');
            fwrite(fid_write,fread(fid_read,3,'int16'),'int16');
        if mod(k,round(n_samples/n_skip_accel/100))==1
            display(['Compass for file ' filename_read ' parsed at ' num2str(k/n_skip_accel/n_samples*100,'%02.2f') '% in ' num2str(toc) 's.']);
        end
        end
end

time_compass_read=toc;
tic;

%% write conductivity data
switch type
    case {}%'v5','v5b'}
        fwrite(fid_write,hex2dec('ts4F'),'uint16');%cond_data tag 'ts4F'
        for k=1:length(pos_cond)
            fseek(fid_read,double(pos_cond(k)),'bof');
            fwrite(fid_write,fread(fid_read,1,'uint32'),'uint32');
        if mod(k,round(n_samples/n_skip_accel/100))==1
            display(['Conductivity for file ' filename_read ' parsed at ' num2str(k/n_skip_accel/n_samples*100,'%02.2f') '% in ' num2str(toc) 's.']);            
        end        
        end
end

time_cond_read=toc;
tic;

%% write message errors
% fwrite(fid_write,hex2dec('F0FF'),'uint16');%msg_data tag 'F0FF'
% for pos=pos_msg
%     fseek(fid_read,pos,'bof');
%     fwrite(fid_write,fread(fid_read,3,'uint32'),'uint32');
% end
% 
% fseek(fid_write,0,'eof');

%% write footer
% fwrite(fid_write,hex2dec('0FFF'),'uint16');%footer tag '0FFF'
% fseek(fid_read,-108,'eof');
% fwrite(fid_write,fread(fid_read,108));
% fseek(fid_write,-8,'eof');
% 
% % p is used instead of n_samples (cheating...)
% fwrite(fid_write,p,'uint32');
% 

fclose(fid_read);

time_total=time_temp_read+time_accel_read+time_compass_read+time_cond_read+toc;

display(['File ' filename_read ' converted in ' num2str(time_total) 's.']);            





