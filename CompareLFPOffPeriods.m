dirname = 'D:\OCD-EEG-LFP-Match\';
%{
filenames = {...
    'aDBS002/2019-06-06/MSIT/aDBS002_MSIT_2019-06-06_14-41-42_synced_eeg_lfp.mat', ...
    'aDBS002/2019-06-06/interview/aDBS002_interview_2019-06-06_11-22-55_synced_eeg_lfp.mat', ...
    'aDBS003/2019-05-16/MSIT/aDBS003_MSIT_2019-05-16_14-04-11_synced_eeg_lfp.mat', ...
    'aDBS003/2019-05-16/beads/aDBS003_beads_2019-05-16_13-41-03_synced_eeg_lfp.mat', ...
    'aDBS003/2019-05-16/interview/aDBS003_interview_2019-05-16_11-12-28_synced_eeg_lfp.mat', ...
    'aDBS003/2019-05-16/programming/aDBS003_programming_2019-05-16_15-30-23_synced_eeg_lfp.mat', ...
    'aDBS003/2019-06-17/MSIT/aDBS003_MSIT_2019-06-17_14-18-05_synced_eeg_lfp.mat', ...
    'aDBS003/2019-06-17/beads/aDBS003_beads_2019-06-17_13-48-38_synced_eeg_lfp.mat', ...
    'aDBS003/2019-06-17/interview/aDBS003_interview_2019-06-17_11-15-19_synced_eeg_lfp.mat', ...
    'aDBS003/2019-06-17/resting-state/aDBS003_resting-state_2019-06-17_13-33-34_synced_eeg_lfp.mat', ...
    'aDBS003/2019-06-17/resting-state/aDBS003_resting-state_2019-06-17_16-22-08_synced_eeg_lfp.mat'};
NumDataSets = numel(filenames);
subplotC = 4;
subplotR = ceil(NumDataSets/subplotC);
for DataLoop = 1:NumDataSets
    fname = filenames{DataLoop};

    fprintf(['\n\nDataSet ' num2str(DataLoop) '\n' fname '\nloading data ... ']); tic

    load([dirname fname])
    LFPlt = lfp_match.left;
    figure
    std_grad = movstd(gradient(LFPlt), 4);
    plot(std_grad)
    cutoff = input("Input STD cutoff: ");
    figure
    hold on
    range = 1:length(LFPlt);
    plot(LFPlt)
    plot(range(std_grad<cutoff), LFPlt(std_grad < cutoff))
end


start = SegLFP0(1)+10000;
stop =  SegLFP0(2)-1000;
%}
fig = uifigure;
ax = uiaxes(fig);
% fig2 = uifigure;
% ax2 = uiaxes(fig2);
max = 1.334010000000000;
min = 1.333990000000000;
sld1 = uislider(fig,...
    'Position',[100 75 120 3],...
    'ValueChangedFcn',@(sld,event) updatePeriod(sld, LFPlt, ax),...
    'Limits', [min,max]);

% sld2 = uislider(fig2,...
%     'Position',[200 75 120 3],...
%     'ValueChangedFcn',@(sld,event) updateBounds(sld, LFPlt, PeriodLFP, ax, ax2),...
%     'Limits', [1, length(LFPlt)]);

% 1.334534000000000
% 1.329348000000000
% 1.331337480000000 G
% 1.332002000000000
% 1.332335000000000
% 1.332477000000000
% 1.332668000000000 M
% 1.332934000000000 G
% 1.333111000000000 G
% 1.333333000000000 M
% 1.333619000000000
% 1.333999960000000 G

%sld2.Value = SegLFP0(2);
sld1.Value = 1.334000000000000;
updatePeriod(sld1, LFPlt, ax)
%{
for i=1:1000
    sld1.Value = min + (max-min)*i/1000;
    updatePeriod(sld1, LFPlt, ax)
    pause(0.1)
end
%}
function updatePeriod2(per, LFPlt)
for i=1:length(pers)
    per = pers(i)
    start = 1*10^5;
    stop = 1.5*10^5;
    inds = start:stop;
    inds = mod(inds, per);
    [inds, ord] = sort(inds);
    temp = LFPlt(start:stop);
    plot(inds, smooth(circshift(temp(ord), floor(length(temp)/2)), 600))
end
end
function updatePeriod(sld, LFPlt, ax)
    start = 1*10^5;
    stop = 1.5*10^5;
    inds = start:stop;
    inds = mod(inds, sld.Value);
    plot(ax, inds, LFPlt(start:stop), 'o')
end
%{
function updateBounds(sld, LFPlt, periodLFP, ax, ax2)
    start = 1*10^5;
    stop = floor(sld.Value);
    samples = start:stop;
    plot(ax, LFPlt)
    hold(ax)
    plot(ax, samples, LFPlt(samples))
    hold(ax)
    xlim(ax, [stop-200, stop+200])
    
    aligned_samples = mod(samples, periodLFP);
    lfp_samples = LFPlt(samples);
    [sorted_samples, I] = sort(aligned_samples);
    sorted_lfp = lfp_samples(I);
    plot(ax2, sorted_samples, sorted_lfp, 'o')
    hold(ax2)
    plot(ax2, sorted_samples, movmin(sorted_lfp, 10))
    plot(ax2, sorted_samples, movmax(sorted_lfp, 10))
    hold(ax2)
end
%}