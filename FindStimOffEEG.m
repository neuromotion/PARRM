function SegEEG = FindStimOffEEG(EventTimesEEG,numEEG,intL)

if nargin < 3 || isempty(intL), intL = 50; end

% intervals between events
EventTimesEEG = [0;EventTimesEEG(:);numEEG+1];
IntervalsEEG = diff(EventTimesEEG);

% find long intervals
EEGlong = find(IntervalsEEG > intL*median(IntervalsEEG));
SegEEG = [EventTimesEEG(EEGlong),EventTimesEEG(EEGlong+1)];

% merge consecutive intervals
pruned = true;
while pruned
    pruned = false;
    for k = 1:size(SegEEG,1)-1
        if SegEEG(k,2) == SegEEG(k+1,1)
            SegEEG(k,2) = SegEEG(k+1,2);
            SegEEG(k+1,:) = [];
            pruned = true;
            break
        end
    end
end
