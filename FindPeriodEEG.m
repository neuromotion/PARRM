function [Period,Seg] = FindPeriodEEG(EEG,EventTimes)

minI = 20;
maxI = 2000;
tolP = .02;
optN = [5000,10000,25000];
midN = [.4;.3;.05];
bw = [50;100;150];

medEEG = diff(EEG); 

IntervalsEEG = diff(EventTimes);

Period = mean(IntervalsEEG(IntervalsEEG > minI & IntervalsEEG < maxI));
regularNdx = (IntervalsEEG > (1-tolP)*Period & IntervalsEEG < (1+tolP)*Period);
Period = mean(IntervalsEEG(regularNdx));
regularNdx = (IntervalsEEG > (1-tolP)*Period & IntervalsEEG < (1+tolP)*Period);
for k = 2:5
    regularNdx = regularNdx | (IntervalsEEG > (k-tolP)*Period & IntervalsEEG < (k+tolP)*Period);
end

%---------------------------------
% optimize the EEG period
%---------------------------------

% find a long segment of regular stimulation
gapIntervals = [0;find(~regularNdx);numel(regularNdx)+1];
[~,k] = max(diff(gapIntervals));
segStart = sum(IntervalsEEG(1:gapIntervals(k)))+1;
segEnd = sum(IntervalsEEG(1:gapIntervals(k+1)));
segLen = segEnd-segStart+1;

Seg = [segStart,segEnd];

% optimize
optN = unique(min(segLen,optN));

for k = 1:numel(optN)
    
    a = segStart+floor(midN(k)*segLen);
    b = segEnd-ceil(midN(k)*segLen);
    t = unique(randperm(b-a,min(optN(k),b-a))+a-1);
    X = [t(:),medEEG(t)];
    Period = fminsearch(@(p)sum(lfpreg(X,p,bw(k)).^2),Period,optimset('display','off'));
    
end
