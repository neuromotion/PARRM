%% FindPeriodLFP.m
% Finds the exact period of a stimulation artifact in an LFP recordinng
%
%% Inputs:
%
%   LFP     : the data in which to find the period of stimulation
%   Seg     : the span in samples to use for finding the period
%   Period  : an initial guess for the value of the period
%% Outputs:
%
%   Period  : an exact estimate for the period of stimulation

function Period = FindPeriodLFP(LFP,Seg,Period)

optN = [5000,10000,25000];
useC = [0,0,.95];

bw = [5,10,20];
lam = 1;

% boundary for first derivative outliers (maxC x mean)
maxC = 3;

%---------------------------------
% standardize
%---------------------------------

if iscell(LFP)
    for k = 1:numel(LFP)
        % work with first derivatives of LFP
        LFP{k} = diff(LFP{k}); 
        % standardize
        LFP{k} = LFP{k}(:) / mean(abs(LFP{k}));
        % clip outliers
        LFP{k} = max(min(LFP{k},maxC),-maxC);
    end
else
    LFP = diff(LFP); 
    LFP = LFP(:) / mean(abs(LFP));
    LFP = max(min(LFP,maxC),-maxC);
end

%---------------------------------
% optimize the LFP period
%---------------------------------

segStart = Seg(1);
segEnd = Seg(2);
segLen = segEnd-segStart+1;

% optimize
optN = unique(min(segLen,optN));

% suppress non-invertible warnings
warning off MATLAB:singularMatrix
warning off MATLAB:nearlySingularMatrix

% Finding the LFP period is tricky because it is very sensitive to slight
% changes in the period and because there are many distractor solutions,
% some very(!) close to the "true" solution, that give a good fit, but that 
% have multiple peaks. This seems to be the result of some weird phenomenon
% in which certain periods mimic certain other periods over long intervals.
%
% This next code block uses a grid search to narrow down to some sensible
% initial values for Matlab's built optimizer. It also begins with smaller
% segments of data to reduce the sensivity of the object to the period and
% progresses to larger segments. It uses fewer basis functions when working
% smaller segments and it also penalizes solutions that need higher- 
% frequency cosines and sines for fitting. This last trick helps avoid
% distractor solutions that find periods that mimic multiples of the "true"
% period.


for k = 1:numel(optN)
    
    % Begin with a stretch of values optN(k) long centered in the segment
    a = ceil((segStart+segEnd-optN(k))/2);
    b = floor((segStart+segEnd+optN(k))/2);
    if segLen*useC(k) < b-a
        t = (a:b).';
    else
        a = segStart+floor((1-useC(k))/2*segLen);
        b = segEnd-ceil((1-useC(k))/2*segLen);
        t = unique(randperm(b-a,min(optN(k),b-a))+a-1);
        t = t(:);
    end

    % generate range of periods to test
    periods = [];
    for j = 1:numel(Period)
        periods = [periods,Period(j)*(1+(-1e-2:1e-4:1e-2)/k),Period(j)*(1+(-1e-3:1e-5:1e-3)/k)];
    end
    periods = unique(periods);
    np = numel(periods);
    v = zeros(np,1);
    for j = 1:np
        v(j) = opt_local(t,LFP,periods(j),min(bw(k),floor(numel(t)/4)),lam);
    end
    [v,vn] = sort(v);
    periods = periods(vn);
    for j = 1:min(5,np)
        [periods(j),v(j)] = fminsearch(@(p)opt_local(t,LFP,p,min(bw(k),floor(numel(t)/4)),lam),periods(j),optimset('display','off'));
    end
    
    [~,j] = min(v);
    Period = periods(j);
end

% do one final optimization with 0 regularization to get the best fit
Period = fminsearch(@(p)opt_local(t,LFP,p,min(bw(k),floor(numel(t)/4)),0),Period,optimset('display','off'));

% turn warnings back on
warning on MATLAB:singularMatrix
warning on MATLAB:nearlySingularMatrix

return

%--------------------------------------------------%
%------------- HELPER FUNCTION --------------------%
%--------------------------------------------------%

function f = opt_local(t,LFP,p,bw,lam)

f = 0;

s = (1:(2*bw+1)); s = lam*s/sum(s);

if iscell(LFP)
    for k = 1:numel(LFP)
        [r,b] = lfpreg([t,LFP{k}(t)],p,bw);
        f = f + mean(r.^2) + s*(b.^2);
    end
    f = f / numel(LFP);
else
    [r,b] = lfpreg([t,LFP(t)],p,bw);
    f = mean(r.^2) + s*(b.^2);
end
