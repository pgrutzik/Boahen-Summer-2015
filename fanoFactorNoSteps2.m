function [FF] = fanoFactorNoSteps2(consideredSpikes, startTime, stopTime) %INSTEAD JUST IMPUT CONSIDERED SPIKES this should be conspikes{n}{y}{c}
%countOnOff is 1 if you want to count the number of spikes in only the on
%and off states separately% spikeWindowLength is in seconds
%idx is the index of conditions you want to analyze
% onOffPhases will be a string. Either 'yes' if you want to calculate the
% Fano factor for each on off phase, or 'no' if you want to calculate the
% Fano factor for the entire channel.

numChann  = size(consideredSpikes,1);
numTrial = size(consideredSpikes,2);

FF = zeros(numChann,1);
 
    for iChann = 1:numChann  % for each channel
        numSpikes = zeros(numTrial,1);
            for iTrial = 1:numTrial %length(consideredSpikes{conditionIdx}) % for each trial                
                spikes = consideredSpikes{iChann, iTrial};
                numSpikes(iTrial) = sum(spikes >=startTime(iChann) & spikes <= stopTime(iChann));               
            end    
            FF(iChann) = var(numSpikes)/mean(numSpikes);  
    end
end
