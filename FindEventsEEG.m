function [EventTimes,UseMin] = FindEventsEEG(EEG)

medC = 50;
maxC = 50;
minT = -500;
maxT = 500;
minI = 20;
maxI = 2000;
tolP = .02;

nEEG = numel(EEG);

% add insignificant randomness to break ties
EEG = EEG + (10^5*rand(size(EEG))+1).*max(eps(EEG),realmin);

% median centering
medEEG = EEG-movmedian(EEG,medC);

%---------------------------------
% find EEG peaks
%---------------------------------

% minimum and maximum filtering
minEEG = movmin(medEEG,maxC);
maxEEG = movmax(medEEG,maxC);

% times of EEG peaks
minTimes = find((medEEG==minEEG) & (medEEG < minT)); minTimes(minTimes == 1 | minTimes == nEEG) = [];
maxTimes = find((medEEG==maxEEG) & (medEEG > maxT)); maxTimes(maxTimes == 1 | maxTimes == nEEG) = [];

% intervals between EEG peaks
minIntervals = diff(minTimes);
maxIntervals = diff(maxTimes);

% get typical interval
Period = mean([minIntervals(minIntervals > minI & minIntervals < maxI); ...
    maxIntervals(maxIntervals > minI & maxIntervals < maxI)]);

% repeat min/max filtering to get exactly one peak in a window
minConst = 2*ceil(.75*Period)+1;

% minimum and maximum filtering
minEEG = movmin(medEEG,minConst);
maxEEG = movmax(medEEG,minConst);

% times of EEG peaks
minTimes = find((medEEG==minEEG) & (medEEG < minT)); minTimes(minTimes == 1 | minTimes == nEEG) = [];
maxTimes = find((medEEG==maxEEG) & (medEEG > maxT)); maxTimes(maxTimes == 1 | maxTimes == nEEG) = [];

% intervals between EEG peaks
minIntervals = diff(minTimes);
maxIntervals = diff(maxTimes);

% histogram of intervals
[~,minCount] = uniquecount(minIntervals);
[~,maxCount] = uniquecount(maxIntervals);

% get typical interval for min
Period = mean(minIntervals(minIntervals > minI & minIntervals < maxI));
minNdx = minIntervals > (1-tolP)*Period & minIntervals < (1+tolP)*Period;
minMAD = mad(minIntervals(minNdx));
for k = 2:5
    minNdx = minNdx | (minIntervals > (k-tolP)*Period & minIntervals < (k+tolP)*Period);
end

% get typical interval for max
Period = mean(maxIntervals(maxIntervals > minI & maxIntervals < maxI));
maxNdx = (maxIntervals > (1-tolP)*Period & maxIntervals < (1+tolP)*Period);
maxMAD = mad(maxIntervals(maxNdx));
for k = 2:5
    maxNdx = maxNdx | (maxIntervals > (k-tolP)*Period & maxIntervals < (k+tolP)*Period);
end

% see if there are a lot more max or min peaks
if sum(minNdx) > 1.01*sum(maxNdx)
    UseMin = true;
elseif sum(maxNdx) > 1.01*sum(minNdx)
    UseMin = false;
elseif sum(minNdx) > 1.001*sum(maxNdx)
    UseMin = true;
elseif sum(maxNdx) > 1.001*sum(minNdx)
    UseMin = false;
elseif maxMAD >= minMAD
    UseMin = true;
else
    UseMin = false;
end

if UseMin
    EventTimes = minTimes;
else
    EventTimes = maxTimes;
end

