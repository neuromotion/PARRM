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
%}
filenames = {'aDBS003/2019-06-17/MSIT/aDBS003_MSIT_2019-06-17_14-18-05_synced_eeg_lfp.mat'};
NumDataSets = numel(filenames);

subplotC = 4;
subplotR = ceil(NumDataSets/subplotC);

doPlots = true;
doLoad = false;

rng('default');

tic

% initial guess at periods for 200,30000 fs
PeriodLFPguess = 1.331337;


if doLoad
    load MattScript_11_05_2019.mat
else
    EventTimesEEG = cell(NumDataSets,1);
    PeriodEEG = zeros(NumDataSets,1);
    UseMinEEG = false(NumDataSets,1);
    SegEEG = zeros(2,NumDataSets);
    OffSegsEEG = cell(NumDataSets,2);
    OffSegsLFP0 = cell(NumDataSets,2);
    LFPOffSet0 = zeros(NumDataSets,1);
    LFPScale0 = zeros(NumDataSets,1);
    SegLFP0 = zeros(2,NumDataSets);
    PeriodLFP = zeros(NumDataSets,1);
    AlignWarning = false(NumDataSets,1);
    LFPScale = zeros(NumDataSets,1);
end


% loop over datasets
for DataLoop = 1:NumDataSets
    
    %--------------------------------------------------------------
    % load data
    %--------------------------------------------------------------
 
    fname = filenames{DataLoop};

    fprintf(['\n\nDataSet ' num2str(DataLoop) '\n' fname '\nloading data ... ']); tic

    load([dirname fname])
    
    fprintf(['(time: ' num2str(round(toc)) 's)\n']);
    
    % copy into more convenient variables
    fsLFP = lfp_match.lfp_fs;
    fsEEG = lfp_match.eeg_fs;
    EEG = lfp_match.eeg(:);
    LFPrt = lfp_match.right(:);
    LFPlt = lfp_match.left(:);
    
    numEEG = numel(EEG);
    numLFP = numel(LFPrt);
    
    clear lfp_match
    
    %----------------------------------------------------------------
    % find EEG events
    %----------------------------------------------------------------
    
    fprintf('finding EEG peaks ... '); tic
    
    if ~doLoad, [EventTimesEEG{DataLoop},UseMinEEG(DataLoop)] = FindEventsEEG(EEG); end
    
    if UseMinEEG(DataLoop), fprintf(' using minimum ... '); else, fprintf(' using maximum ... '); end
    fprintf(['(time: ' num2str(round(toc)) 's)\n']);
    
    %----------------------------------------------------------------
    % find EEG period
    %----------------------------------------------------------------
    
    fprintf('optimizing EEG period ... '); tic
    
    if ~doLoad, [PeriodEEG(DataLoop),SegEEG(:,DataLoop)] = FindPeriodEEG(EEG,EventTimesEEG{DataLoop}); end
    
    fprintf(['period=' num2str(PeriodEEG(DataLoop),16) ' ... '])
    fprintf(['(time: ' num2str(round(toc)) 's)\n']);
    
    %----------------------------------------------------------------
    % find approximate intervals of stimulation off
    %----------------------------------------------------------------
    
    fprintf('finding stim-off (approx) ... '); tic
        
    a = ceil(PeriodEEG(DataLoop)*200);
    b = ceil(a*fsLFP/fsEEG);
    
    if ~doLoad
        OffSegsEEG{DataLoop,1} = FindStimOffEEG(EventTimesEEG{DataLoop},numEEG);
        OffSegsEEG{DataLoop,2} = FindStimOff(EEG,a);
        OffSegsLFP0{DataLoop,1} = FindStimOff(LFPrt,b);
        OffSegsLFP0{DataLoop,2} = FindStimOff(LFPlt,b);
    end

    c = [size(OffSegsEEG{DataLoop,1},1),size(OffSegsEEG{DataLoop,2},1),size(OffSegsLFP0{DataLoop,1},1),size(OffSegsLFP0{DataLoop,2},1)];
    fprintf(['# of stim off intervals=(' num2str(c(1)) ',' num2str(c(2)) ',' num2str(c(3)) ',' num2str(c(4)) ') ... '])
    fprintf(['(time: ' num2str(round(toc)) 's)\n']);

    %----------------------------------------------------------------
    % find a rough alignment using stimulation off intervals
    %----------------------------------------------------------------
   
    fprintf('aligning stim-off (approx) ... '); tic
   
    if ~doLoad
        [LFPOffSet0(DataLoop),LFPScale0(DataLoop)] = AlignStimOff(OffSegsEEG(DataLoop,:),OffSegsLFP0(DataLoop,:),fsEEG,fsLFP);
    end
    fprintf(['approx LFP scale: ' num2str(LFPScale0(DataLoop)) ' ... '])
    fprintf(['(time: ' num2str(round(toc)) 's)\n']);

    %----------------------------------------------------------------
    % find LFP period
    %----------------------------------------------------------------

    fprintf('optimizing LFP period ... '); tic
    
    % find the corresponding LFP segment for period determination
    
    SegLFP0(:,DataLoop) = round((SegEEG(:,DataLoop)-LFPOffSet0(DataLoop))/LFPScale0(DataLoop));
    if SegLFP0(1,DataLoop) < 0 || SegLFP0(2,DataLoop) > numLFP+1
        fprintf(['WARNING: Possible misalignment (' num2str(SegLFP0(1,DataLoop)/numLFP,3) ',' num2str(SegLFP0(2,DataLoop)/numLFP,3) ') ... '])
        AlignWarning(DataLoop) = true;
    end
    SegLFP0(:,DataLoop) = max(min(SegLFP0(:,DataLoop),numLFP),1);

    % optimize period
    
    if ~doLoad 
        periods0 = [PeriodEEG(DataLoop)./LFPScale0(DataLoop),200*fsLFP/fsEEG,PeriodLFPguess];
        PeriodLFP(DataLoop) = FindPeriodLFP({LFPlt,LFPrt},SegLFP0(:,DataLoop),periods0);
    end
    
    fprintf(['period=' num2str(PeriodLFP(DataLoop),16) ' ... '])
    fprintf(['(time: ' num2str(round(toc)) 's)\n']);
    
    %----------------------------------------------------------------
    % find the corresponding LFP segment for period determination
    %----------------------------------------------------------------

    fprintf('finding LFP scaling ... '); tic

    scales = PeriodEEG(DataLoop)/PeriodLFP(DataLoop)*[1./[2:10],1:10];
    [~,k] = min(abs(scales-LFPScale0(DataLoop)));
    LFPScale(DataLoop) = scales(k);
    
    fprintf(['LFP Scale=' num2str(LFPScale(DataLoop),16) ' ... '])
    fprintf(['(time: ' num2str(round(toc)) 's)\n']);
    
    % plotting
    if doPlots
        
        figure(1), set(gcf,'name','EEG')
        
        subplot(subplotR,subplotC,DataLoop)
        
        segStart = SegEEG(1,DataLoop); segEnd = SegEEG(2,DataLoop); segLen = segEnd-segStart+1;
        
        medEEG = EEG-movmedian(EEG,50);
        
        t = unique(randperm(segLen,min(segLen,50000))+segStart-1);
        X = [t(:),medEEG(t)];
        [~,~,f,fx] = lfpreg(X,PeriodEEG(DataLoop),150);
        
        plot(mod(t,PeriodEEG(DataLoop)),medEEG(t),'k.',fx,f,'y-','linewidth',1)
        title(['DataSet ' num2str(DataLoop) ' (' num2str(round(segLen/10^6)) ' million bins)'])
        set(gca,'xlim',[0 PeriodEEG(DataLoop)])
        
        figure(2), set(gcf,'name','LFP')
        
        subplot(subplotR,subplotC,DataLoop)
        
        segStart = SegLFP0(1,DataLoop); segEnd = SegLFP0(2,DataLoop); segLen = segEnd-segStart+1;
        
        medLFP = LFPlt-movmedian(LFPlt,667);
        medLFP = medLFP / mean(abs(medLFP));
        
        t = segStart+floor(.2*segLen):segEnd-ceil(.2*segLen);
        X = [t(:),medLFP(t)];
        [~,~,f,fx] = lfpreg(X,PeriodLFP(DataLoop),150);
        
        plot(mod(t,PeriodLFP(DataLoop)),medLFP(t),'b.',fx,f,'y-','linewidth',1)
        hold on
        
        medLFP = LFPrt-movmedian(LFPrt,667);
        medLFP = medLFP / mean(abs(medLFP));
        
        X = [t(:),medLFP(t)];
        [~,~,f,fx] = lfpreg(X,PeriodLFP(DataLoop),150);
        
        plot(mod(t,PeriodLFP(DataLoop))+PeriodLFP(DataLoop),medLFP(t),'r.',fx+PeriodLFP(DataLoop),f,'g-','linewidth',1)
        hold off

        title(['DataSet ' num2str(DataLoop) ' (' num2str(round(segLen/1000)) ' thousand bins)'])
        set(gca,'xlim',[0 2*PeriodLFP(DataLoop)],'ylim',[-5 5])
        
        drawnow, pause(.1)
    end
    
end

if false
    save MattScript_11_05_2019.mat dirname filenames EventTimesEEG PeriodEEG ...
        UseMinEEG SegEEG OffSegsEEG OffSegsLFP0 LFPOffSet0 LFPScale0 SegLFP0 ...
        PeriodLFP AlignWarning LFPScale
end


% % adjust approximate intervals
% for k = 1:size(OffSegsEEG,1)
%     a = OffSegsEEG(k,1);
%     b = OffSegsEEG(k,2);
%     
%     overlap = find( ((SegEEG(:,1) <= a) & (a <= SegEEG(:,2))) ...
%         | ((SegEEG(:,1) <= b) & (b <= SegEEG(:,2))) ...
%         | ((a <= SegEEG(:,1)) & (SegEEG(:,1) <= b)) ...
%         | ((a <= SegEEG(:,2)) & (SegEEG(:,2) <= b)) );
%     
%     if ~isempty(overlap)
%         OffSegsEEG(k,1) = min(SegEEG(overlap,1));
%         OffSegsEEG(k,2) = max(SegEEG(overlap,2));
%     end
% end