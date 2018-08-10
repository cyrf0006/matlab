function quickLookFSI(file, varargin)
    
% function quickLookFSI(file)
%
% This function uses modified files from the .xlsx I recieved from
% Furu. I basicaly copy-pasted the content without the header and I
% replaced '/' and ':' by spaces in a shell or Emacs.
%
% Columns should be as follow:
%   1:5 - time without second (useless)
%     6 - AVN (north vel.)
%     7 - AVE (east vel.)
%     8 - AVU (vert. vel)
%     9 - ASPD (average speed)
%    10 - AVDIR (direction 0-360deg. rel. to north)
%    11 - ATLT (intrument tilt in degree)
% 12:17 - TIME (ss mm hh MM DD YYY)
%    18 - DIR
%    19 - BATT
%    20 - TX
%    21 - TY
%    22 - VN
%    23 - VE
%    24 - VU
%    25 - STEMP
%    26 - COS(DIR) (if present...)   
%
% usage ex: 
% >> quickLookFSI('FSI_041.dat', [datenum(2012,10,11,0,0,0) datenum(2012,10,13,6,30,0)])
% >> quickLookFSI('FSI_005.dat', [datenum(2012,10,7,14,0,0) datenum(2012,10,10,7,0,0)])  
    

    
data = load(file);

timeVec = datenum(data(:,17), data(:,15), data(:,16), data(:,12), ...
                  data(:,13), data(:,14));

if ~isempty(varargin)
    range = varargin{:};
    if length(range) == 2
        I = find(timeVec>=range(1) & timeVec<=range(2));
        timeVec = timeVec(I);
        data = data(I,:);
    else
        disp('problem with varargin')
        return
    end
end

E = data(:,7);
N = data(:,6);
W = data(:,8);
T = data(:,25);

[U,V] = rotate_vecd(E,N,-30);


figure(1)
clf
plot(timeVec, V)
datetick
hold on
plot(timeVec, U, 'r')

figure(2)
clf
plot(timeVec, T)
datetick



keyboard