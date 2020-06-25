# PARRM

These are all the files we've used over the course of developing the method. I've included them all here in case there's interest but not all of them are necessary for the current incarnation. 

This method operates on the assumption that the artifact is semi-regular, periodic, and linearly combined with the signal of interest.
The first step is finding the period of the artifact in the raw signal. This is done using a grid search centered about the quotient of the sampling rate and stimulation frequency. In order to evaluate each period, samples are translated onto the timescale of the period by finding the modulo of each sample number with the current period of interest. Then, using linear regression, the best fitting sum of sinusoidal harmonics of the period is found and the mean squared error of the fit with the raw data is computed. The period which minimizes this mean squared error is then chosen. This period is used to produce a linear filter that takes four parameters in addition to the period itself: the window size, the skip size, the period distance, and the window direction. Only samples within the window size are used to filter the data at each sample. Samples within the skip size are not used for filtering. Only samples within the period distance on the timescale of the period are averaged together. The window direction is used to specify whether the window only takes past, future, or both into account.

This is typically done using the following lines of code:

'''Matlab
% Filter timeseries 'data' using PARRM
% Assume 'data' has a 200Hz sampling rate and 150Hz stimulation frequency
guessPeriod=200/150; % Theoretical stimulation period
span=[2000,12000]; % Span of samples in 'data' where artifact is regular
windowSize=2000; % Width of window in sample space to be used for removal
skipSize=20; % Number of samples to ignore in sample space
windowDirection='both'; % Remove using information from the past and future

Period=FindPeriodLFP(data,span,guessPeriod); % Find the period of stimulation in 'data'
periodDist=Period/120; % Window in period space for which samples will be averaged

PARRM=PeriodicFilter(Period,windowSize,periodDist,skipSize,windowDirection); % Create the linear filter
Filtered=((filter2(PARRM.',data,'same')-data)./(1-filter2(PARRM.',ones(size(data)),'same'))+data)'; % Filter using the linea filter and remove edge effects
'''
