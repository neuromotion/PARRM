function s = AlignScore(x,LFPSegs,EEGSegs)

LFPSegs = x(1)+x(2)*LFPSegs;

s = 0;
for k = 1:size(LFPSegs,1)
    s = s + sum(max(min(LFPSegs(k,2),EEGSegs(:,2))-max(LFPSegs(k,1),EEGSegs(:,1)),0)); 
end
s = 0.5*(s/(size(EEGSegs,1)*sum(LFPSegs(:,2)-LFPSegs(:,1)))+s/(size(LFPSegs,1)*sum(EEGSegs(:,2)-EEGSegs(:,1))));












