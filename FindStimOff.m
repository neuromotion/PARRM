function Segs = FindStimOff(Data,WindowWidth,MinWidth,FillWidth)
% function Segs = FindStimOff(Data,WindowWidth,MinWidth,FillWidth)

if nargin < 3 || isempty(MinWidth)
    MinWidth = WindowWidth;
end
if nargin < 4 || isempty(FillWidth)
    FillWidth = ceil(MinWidth/2);
end

% get fluctuations around median
Data = abs(Data - movmedian(Data,ceil(WindowWidth/4)));

% get maximum fluctuations in a region near each time point
% and accentuate regions with small maximum fluctuations
v = log(min(movmax(Data,[WindowWidth 1]),movmax(Data,[1 WindowWidth]))+1);
v = min(v,median(v));

% find points with low fluctuations
gap = v < max(v)-4*std(v);
N = numel(gap);

% find regions that tend to have low fluctuations
on = false;
j = 0;
for k = 1:N
    if gap(k)
        if on
            gap(j:k) = true;
        else
            on = true;
        end
        j = k;
    elseif on && k-j > FillWidth
        on = false;
    end
end

% prune regions that are too small
on = false;
j = 0;
for k = 1:N
    if gap(k) && ~on
        on = true;
        j = k;
    elseif ~gap(k) && on
        if k-j < MinWidth
            gap(j:k) = false;
        end
        on = false;
    end
end
if on && N+1-j < MinWidth, gap(j:N) = false; end

% get the starting and ending of segments with low fluctations
D = diff([false;gap;false]);
Segs = [find(D==1),find(D==-1)-1];




