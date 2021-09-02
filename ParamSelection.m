%% Import data and find period
load('simData')
raw=simData.downedRecording; % Vector containing neural data
fs=200; % Sampling rate
stimRate=150; % Stimulation frequency
Period=FindPeriodLFP(raw,[1,length(raw)-1],fs/stimRate); % Find period
%% Param tester
% Plots the result of a PARRM filter for a window size and period distance
% All samples are plotted on the timescale of the period on the left
% Raw and filtered timeseries are plotted on the right
% Zoom in horizontally to select your period distance
% Use the slider in the bottom left to change your window size

initWinSize=2000; % Initialize value for the window size parameter
figure('units','normalized','outerposition',[0 0 1 1])
ax1=subplot(1,2,1);
plot(mod(1:initWinSize*2-1,Period),diff(raw(1:initWinSize*2)),'o')
axis tight
ax2=subplot(1,2,2);
plot(raw)
axis tight
h2 = uicontrol('style','slider','units','pixel','position',[20 20 300 20],'min',1,'max',round((length(raw)-1)/2),...
'callback',@(hObject, event) reset(ax1,ax2,raw,Period,hObject),'val',initWinSize);
h=zoom;
h.ActionPostCallback = @(obj,evd) zoomcallback(obj,raw,Period,h2);
h.Enable = 'on';
h3=uicontrol('style','text','position',[20 40 300 20],'String',"Window size ");
reset(ax1,ax2,raw,Period,h2)

function zoomcallback(obj,raw,Period,sld)
    allAxes = findobj(obj.Children,'Type','Axes');
    numClicked = find(gca==allAxes);
    if numClicked==2
        reset(allAxes(2),allAxes(1),raw,Period,sld)
    end
end

function reset(ax1,ax2,raw,Period,sld)
    newLim = xlim(ax1);
    xl2=xlim(ax2);
    yl2=ylim(ax2);
    hold(ax2,'off')
    plot(ax2,raw)
    PARRM=PeriodicFilter(Period,ceil(get(sld,'Value')),diff(newLim)/2,0,'both');
    Filtered=((filter2(PARRM.',raw','same')-raw')./(1-filter2(PARRM.',ones(size(raw')),'same'))+raw')';
    hold(ax2,'on')
    plot(ax2,Filtered,'LineWidth',1)
    legend(ax2,["Raw","Filtered"])
    title(ax2,'Data')
    xlabel(ax2,'Sample #')
    xlim(ax2,xl2)
    ylim(ax2,yl2)
    xl=xlim(ax1);
    yl=ylim(ax1);
    plot(ax1,mod(1:2*ceil(get(sld,'Value'))-1,Period),diff(raw(1:2*ceil(get(sld,'Value')))),'o')
    xlim(ax1,xl);
    ylim(ax1,yl);
    xlabel(ax1,'Sample in Period')
    ylabel(ax1,'Amplitude')
    title(ax1,{strcat("WindowSize: ",num2str(ceil(get(sld,'Value')))),strcat("PeriodWidth: ",num2str(diff(newLim)/2))})
end