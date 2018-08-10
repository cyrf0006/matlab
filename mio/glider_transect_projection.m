function  [theIndex, distVec] = glider_transect_projection(lonVec, latVec, origin, target, varargin)

% usage ex:
% xVec = transect_projection(lonVec, latVec, origin, target, timeLims)

if ~isempty(varargin)
    timeVec = varargin{1};
    timeLims = varargin{2};
    I = find(timeVec>=timeLims(1) & timeVec<=timeLims(2));
    lonVec_orig = lonVec;
    latVec_orig = latVec;
    index_orig = 1:length(lonVec);
    lonVec = lonVec(I);
    latVec = latVec(I);
    timeVec = timeVec(I);
    index_orig = index_orig(I);
end

transect_distance = m_lldist([origin(2) target(2)], [origin(1) target(1)]);

[Y, I2] = min(abs(latVec-target(1)) + abs(lonVec-target(2)));
%[Y, I1] = min(abs(latVec(1:I2)-origin(1)) + abs(lonVec(1:I2)-origin(2)));
%theIndex = I1:I2;
[Y, I1] = min(abs(latVec-origin(1)) + abs(lonVec-origin(2)));
theIndex = min([I1 I2]):max([I1 I2]);
[distance, lon, lat] = m_lldist([lonVec(I1) lonVec(I2)], [latVec([I1]) latVec([I2])], 10000);
distVec = nan(size(theIndex));
for i = 1:length(theIndex)
    [Y, I] = min(abs(latVec(theIndex(i))-lat) + abs(lonVec(theIndex(i))-lon));
    distVec(i) = m_lldist([origin(2) lon(I)], [origin(1) lat(I)]);
end

%% Clean trajectory (remove repetition + backward trajectory)
% ---  Version 1  --- %
theIndexUnique = theIndex;
Inegative = find(diff(distVec)<0);
p = polyfit(1:length(distVec), distVec, 1);
if p(1)>0 % origin -> target
    while ~isempty(Inegative)
        Inegative = find(diff(distVec)<0);
        if ~isempty(Inegative)
            distVec(Inegative+1) = [];
            theIndexUnique(Inegative+1) = [];
        end
        Inegative = find(diff(distVec)<0);
    end
else % target -> origin (return trajectory)
    while ~isempty(Inegative)
        Inegative = find(diff(distVec)>0);
        if ~isempty(Inegative)
            distVec(Inegative+1) = [];
            theIndexUnique(Inegative+1) = [];
        end
        Inegative = find(diff(distVec)>0);
    end
end

% $$$ % ---  Version 2  --- %
% $$$ theIndexUnique = theIndex;
% $$$ [distVec, Isort] = sort(distVec);
% $$$ timeVec = timeVec(Isort)

% $$$ distVec(Inegative+1) = [];
% $$$ theIndexUnique(Inegative+1) = [];
%timeVec = timeVec(I);

[distVec_unique, I] = unique(distVec, 'last');
distVec = interp1(theIndexUnique(I), distVec_unique, theIndex);

if ~isempty(varargin)
    theIndex = index_orig(theIndex);
end
