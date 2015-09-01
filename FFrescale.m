function [] = FFrescale(startTime, stopTime)

load('~/Documents/BoahenLabProject15/poolAllData_numState_2.mat');
load('~/Documents/BoahenLabProject15/hash.mat');

[stopStartTimes] = findOnOffTimes(popData); 

recIdx = strcmp(popParam.recording, '2011_12_17');

fieldNames = {'covert', 'overt', 'control'};

uniqueTheta = unique(popParam.theta);

binSize = .01;

for i = 1:length(fieldNames)
    FForig.(fieldNames{i}) = zeros(length(uniqueTheta), size(CueAlign ,1));
    FFresc.(fieldNames{i}) = zeros(length(uniqueTheta), size(CueAlign ,1));
end;

for iTheta = 1:length(uniqueTheta)
    thisTheta = uniqueTheta(iTheta);
    thetaIdx = popParam.theta == thisTheta;
    
    for iField = 1:length(fieldNames)
        switch fieldNames{iField}
            case 'covert'
                condIdx = popParam.iCond == 3;
            case 'overt'
                condIdx = popParam.iCond == 1;
            case 'control'
                condIdx = popParam.iCond == 2 | popParam.iCond == 4;
            otherwise
                error('Write some error');
        end;

        idx = find(thetaIdx & recIdx & condIdx);     
        
        allAllTrial = [];
        consideredRescSpikes = [];
        
        for i = length(idx)      
            allTrial = popParam.behaviorParam{idx(i)}.trialNumbers;
            firstTimeOrig = cellfun( @(x) x(1), popData.timeBins{idx(i)}); % beginning of analyzed spikes
            endTimeOrig = cellfun( @(x) x(end), popData.timeBins{idx(i)});         
            allTrial(endTimeOrig < stopTime) = [];  %%%%%
            
            firingRatesOn = popData.estEmis{idx(i)}(2,:) / binSize; % firing rate in Hz. estEmis is the emission probablility estEmis = firingrate * 0.01
            firingRatesOff = popData.estEmis{idx(i)}(1,:) / binSize;
            
            rescaledSpikes = rescaledSpikeValues2(CueAlign(:, allTrial), stopStartTimes{idx(i)}, firingRatesOn, firingRatesOff);
            
            consideredRescSpikes = [consideredRescSpikes, rescaledSpikes];
            allAllTrial = [allAllTrial, allTrial];
        end  

        firstTimeOrig(1:size(popData.timeBins{idx(i)})) = startTime;
        
        consideredSpikes = CueAlign(:, allAllTrial);   % If unattended, considered spikes is no longer associated with a specific condition number. So, we have to calculate firing rate and re           
        avgFiringRates = getFiringRate2(consideredSpikes, firstTimeOrig, endTimeOrig); % To rescale the time window later
        
        firstTimeResc = zeros(size(popData.timeBins{idx(i)},1));
        avgFiringRatesResc = getFiringRate2(consideredRescSpikes, firstTimeResc, firstTimeResc);
        
        origStartTimes(1:size(consideredSpikes,1)) = startTime; 
        origStopTimes(1:size(consideredSpikes,1)) = stopTime;
        timeWindow = stopTime - startTime;
        rescaleUnits = arrayfun(@(x) x * timeWindow, avgFiringRates);
       
        arrayNumiChann = 1:size(consideredSpikes,1); 
        rescTimeWindow = arrayfun(@(x) rescaleUnits(x)/avgFiringRatesResc(x), arrayNumiChann);        
        rescStartTimes(1:size(consideredRescSpikes,1)) = 0; % Change
        
        FForig.(fieldNames{iField})(iTheta,:) = fanoFactorNoSteps2(consideredSpikes, origStartTimes, origStopTimes);         
        FFresc.(fieldNames{iField})(iTheta,:) = fanoFactorNoSteps2(consideredRescSpikes, rescStartTimes, rescTimeWindow); 
    end;
end

figure;
hold on;
set(plot( mean(( FForig.overt-FForig.control)/(FForig.overt+FForig.control)), 'bo' ),'linewidth',3);
set(gca,'linewidth',3);
x = 0:1:length(mean(( FForig.overt-FForig.control)/(FForig.overt+FForig.control)));
y(1:length(x)) = mean( mean(( FForig.overt-FForig.control)/(FForig.overt+FForig.control)));
set(plot(x,y, 'b--'),'linewidth',3);
set(gca,'linewidth',3);
%plot(mean(mean(( FForig.overt-FForig.control)/(FForig.overt+FForig.control))),'b--')

%xlabel('Unique theta');
%ylabel('Fano Factor');
%title('Overt Modulation Index');
%print('overtFForig', '-dpdf')

set(plot( mean(( FFresc.overt-FFresc.control)/(FFresc.overt+FFresc.control)), 'ro' ),'linewidth',3);
set(gca,'linewidth',3);
x = 0:1:length( mean(( FFresc.overt-FFresc.control)/(FFresc.overt+FFresc.control)));
y(1:length(x)) = mean( mean(( FFresc.overt-FFresc.control)/(FFresc.overt+FFresc.control)))
set(plot(x,y, 'r--'),'linewidth',3);
set(gca,'linewidth',3);

xlabel('Unique theta');
set(gca,'linewidth',3);
ylabel('Fano Factor');
set(gca,'linewidth',3);
title('Overt Fano Factor Modulation Index');
set(gca,'linewidth',3);
legend('Origional Modulation Index','Rescaled ISI Modulation Index')
set(gca,'linewidth',3);
print('overtFFdiff', '-dpdf')

figure;
hold on
set(plot( mean( (FForig.covert-FForig.control)/(FForig.covert+FForig.control)  ), 'ro' ),'linewidth',3);
x = 0:1:length(mean( (FForig.covert-FForig.control)/(FForig.covert+FForig.control)));
y(1:length(x)) = mean( mean( (FForig.covert-FForig.control)/(FForig.covert+FForig.control)  ))
set(plot(x,y, 'b--'),'linewidth',3);
set(gca,'linewidth',3);

% xlabel('Unique theta');
% ylabel('Fano Factor');
% title('Covert Modulation Index');
%print('covertFForig', '-dpdf')

% set(plot( mean( (FFresc.covert-FFresc.control)/(FFresc.covert+FFresc.control)  ), 'bo' ),'linewidth',3);
% set(gca,'linewidth',3);
% x = 0:1:length(mean( (FFresc.covert-FFresc.control)/(FFresc.covert+FFresc.control)  ));
% y(1:length(x)) = mean(mean( (FFresc.covert-FFresc.control)/(FFresc.covert+FFresc.control)  ));
% set(plot(x,y,'r--'),'linewidth',3);
% set(gca,'linewidth',3);

set(plot(mean( FForig.overt), mean(FFmodel.overt),'bo'),'linewidth',3);
set(plot(mean( FForig.covert), mean(FFmodel.covert),'ro'),'linewidth',3);
set(plot(mean( FForig.control), mean(FFmodel.control),'ko'),'linewidth',3);

xlabel('Unique theta');
set(gca,'linewidth',3);
ylabel('Fano Factor');
set(gca,'linewidth',3);
title('Covert Fano Factor Modulation Index');
set(gca,'linewidth',3);
legend('Origional Modulation Index','Rescaled ISI Modulation Index')
set(gca,'linewidth',3);
print('covertFFdiff', '-dpdf')
end    