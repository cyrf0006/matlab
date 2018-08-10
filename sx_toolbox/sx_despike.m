%% despike
% Remove short-duration spikes from a signal.
%%
% <latex>\index{Type A!despike}</latex>
%
%%% Syntax
%   [y, spike] = despike( dv, thresh, smooth, Fs, N )
%
% * [dv] Signal to be despiked.
% * [thresh] Threshold value for ratio of instantaneous signal to rectified, 
%        smoothed signal. A value of 7 is a good starting value. (default: 5)
% * [smooth] Inverse smoothing time-scale (Hz) for the low-pass filtered that 
%        is applied to the high-pass and rectified version of dv to derive an 
%        estimate of the standard deviation of dv. (default: 0.5)
% * [Fs] Sampling frequency (Hz).  (default: 512)
% * [N] Half-width of a spike. The amount of data that gets replaced is 2N+1, 
%        with half on each side of the spike. The replaced data equals the local
%        mean which is estimated after exclusion of the identified spikes. 
%        (default: 5)
% * []
% * [y] Despiked signal.
% * [spike] Index to the identified spikes. It's length is a useful indicator of 
%        the function's veracity.  (optional)
%
%%% Description
% Designed to despike all kinds of signals. It identifies spikes by comparing the
% instantaneous rectified signal against its local standard deviation.  It
% obtains a measure that approximates the local standard deviation by:
% 
% # high-pass filtering the input signal, 
% # rectifying it (absolute value), and
% # smoothing it with a low-pass zero-phase filter.
%
% Based on a program originally written for shear probe despiking (enhance.m).
%
%%% Examples
%
%    >> [y, spike] = despike(dv)
%    >> [y, spike] = despike(dv, thresh)
%    >> [y, spike] = despike(dv, thresh, smooth)
%    >> [y, spike] = despike(dv, thresh, smooth, Fs)
%    >> [y, spike] = despike(dv, thresh, smooth, Fs, N)
%
% @image @images/despike @Despike function @applied with despike(sh1, 7, 0.1, 512, 
% 30). The smoothing scale, the inverse of 0.1 Hz, is long because this profiler 
% was moving at only 0.17 m/s.
% 
% @image @images/despike2 @Despike with zoomed in plot. @The same plot as the 
% previous figure but with a close-up view.
%

% *Version History:*
%
% * 1995-07-27 (FW) initial 
% * 1995-08-22 (FW) revised 
% * 1995-08-28 (FW) revised 
% * 1997-07-14 (FW) revised 
% * 1998-08-31 (RGL) revised 
% * 2000-05-03 (RGL) added unipulse feature 
% * 2000-06-07 (RGL) cured problem with spikes near ends. Improved memory efficiency
% * 2000-06-21 (FW) finally cured problem with spikes near ends.
% * 2002-08-29 (IG) modified spike to correspond to correct location of spikes in y
% * 2009-03-25 (RGL) Set default rate to 512 Hz.
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-05-07 (WID) update to documentation + removed outdated matlab
%                    functions
% * 2012-11-05 (WID) update to documentation
% * 2013-02-26 (WID) merge in changes from Rolf wrt input parameter parsing

function [y, spike] = sx_despike (dv, thresh, smooth, Fs, N)

% check for input parameters
if (nargin > 5), error('Too many input arguments to despike'), end
if (nargin < 1), error('Too few input arguments to despike'), end
if nargin<5, N=6;          end
if nargin<4, Fs = 512;     end
if nargin<3, smooth = 0.5; end
if nargin<2, thresh = 5;   end

dv = dv(:); % force into a column vector

% Zero padding alleviates problems when spikes are at the
% beginning or the end of the vector, and avoids problems
% caused by filter transients of the smoothing filter
% and rectifying filters. (Fab)
len = length(dv);
padLen = 2*fix(Fs/smooth); 
padLen = min([len padLen]);    % padLen must not exceed the length of the input vector
range = [1+padLen len+padLen]; % these two values mark the start and end of the input 
                               % vector in the padded vector (used later)

dv = [flipud(dv(1:padLen)); dv; flipud(dv(len-padLen+1:len))];

% copy vector
y = dv;

% high pass filter & rectify
[b,a] = butter(1, smooth/(Fs/2), 'high');
dv = abs(filtfilt(b,a, dv));

% smoothing coefficients
[b,a] = butter(2, smooth/(Fs/2));

% find spikes
s = warning; % cache warning state 
warning off  % avoid divide-by-zero warnings (fab)
spike = find( dv./filtfilt(b,a,dv) > thresh);
warning(s);  % reset warning state

% ignore spikes detected in the padding
spike(find(spike<range(1) | spike>range(2))) = [];  

if isempty(spike) % get out if no spikes are found (RGL)
   y = y(range(1):range(2));   % remove the padding
   return;
end; 

% exorcise the spikes
% First find spikes that are more than N points apart
if numel(spike) >1
    index = 2;
    while (index <= numel(spike))
        while ((spike(index) - spike(index-1) < N+1) & (numel(spike)>index))
            spike(index) = [];
        end
        index = index + 1;
    end
end

        
        
% On the first pass we will replace the spikes and a range of +-N points around a spike
% with the mean of the adjacent points. This makes sure that the neighborhood mean
% (derived from a low-pass filter)is not biased by a unipolar pulse, such as found in
% conductivity records when bugs go through the sampling volume.
% Unipolar pulses seldom occur in shear probe records.
for k = 1:length(spike); % replace points with mean of points on the edge.
% y(spike(k)-ceil(N/2):spike(k)+N) = (y(spike(k)-ceil(N/2)) + y(spike(k)+N))*ones(1,1+ceil(N/2)+N)/2; 
   select = spike(k)-ceil(N/2):spike(k)+N;
   replace = (y(select(1)) + y(select(end)))/2;
   y(select) = replace;
end

% Now apply low-pass filter with same cut off as specified by input arguments to find a local
% neighborhood mean around the spikes free of bias from unipolar pulses.

% dv = filtfilt(b,a, y); %Filter semi-cleaned version for use in replacement. 
%                        %Won't be biased by large spikes.
% 
% for k = 1:length(spike)
%    select = spike(k)-ceil(N/2):spike(k)+N;
%    y(select) = dv(select);
% %   y(spike(k)-ceil(N/2):spike(k)+N) = dv(spike(k))*ones(1,1+ceil(N/2)+N);
% end

% remove the padding
y = y(range(1):range(2));
spike = spike - padLen;
