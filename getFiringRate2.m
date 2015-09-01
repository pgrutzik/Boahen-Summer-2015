function [avgFiringRate] = getFiringRate2(consideredSpikes, firstTime, endTime)

numChann  = size(consideredSpikes,1);
numTrial = size(consideredSpikes,2);

firingRate = zeros(numChann,numTrial); 
avgFiringRate = zeros(numChann,1);
    for iChann = 1:numChann  % for each channel
            for iTrial = 1:numTrial %length(consideredSpikes{conditionIdx}) % for each trial
                if firstTime(iTrial) == endTime(iTrial)
                    endTime(iTrial) = consideredSpikes{iChann, iTrial}(end);
                end
                spikes = consideredSpikes{iChann, iTrial};
                numSpikes = sum(spikes >= firstTime(iTrial) & spikes <= endTime(iTrial));                
                firingRate(iChann, iTrial) = numSpikes/(endTime(iTrial) - firstTime(iTrial)); 
            end    
            avgFiringRate(iChann) = mean(firingRate(iChann,:));            
    end
end
