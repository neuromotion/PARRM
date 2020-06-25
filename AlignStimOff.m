function [LFPOffSet,LFPScale,fOld] = AlignStimOff(OffSegsEEG,OffSegsLFP,fsEEG,fsLFP)%,numEEG,numLFP)

SegEOff = zeros(0,2);
for k = 1:numel(OffSegsEEG)
    SegEOff = [SegEOff;OffSegsEEG{k}];
end
NE = size(SegEOff,1);

SegLOff = zeros(0,2);
for k = 1:numel(OffSegsLFP)
    SegLOff = [SegLOff;OffSegsLFP{k}];
end
NL = size(SegLOff,2);

c = fsEEG/fsLFP;

if NE*NL == 0
    LFPOffSet = 0; LFPScale = c; fOld = 0;
    return
end

midE = mean(SegEOff,2);
midL = c*mean(SegLOff,2);

h = 0.1;

[xOld,fOld] = fminsearch(@(x)-AlignScore(x,SegLOff,SegEOff)+h*abs(x(2)-c)^2,[0;c],optimset('display','off'));

for k = 1:NE
    for j = 1:NL
        [xNew,fNew] = fminsearch(@(x)-AlignScore(x,SegLOff,SegEOff)+h*abs(x(2)-c)^2,[midE(k)-midL(j);c],optimset('display','off'));
        if fNew < fOld
            fOld = fNew; xOld = xNew;
        end
    end
end

LFPOffSet = xOld(1); LFPScale = xOld(2);




