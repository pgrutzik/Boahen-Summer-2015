function [spikeYvalues] = rescaledSpikeValues2(consideredSpikes, stopStartTimes, firingRatesOn, firingRatesOff)

numChann = size(consideredSpikes, 1);
numTrial = size(consideredSpikes, 2);
spikeYvalues = cell(numChann,numTrial); 
   
    for iTrial = 1:numTrial % for each trial
        
        for iChann = 1:numChann  % for each channel
            
            tOnOffTimes = stopStartTimes{iTrial}(:,1);
            tOnOffBlocks = stopStartTimes{iTrial}(:,2);
            
            start = tOnOffTimes(1);
            stop = tOnOffTimes(end);
            
            numTblock = length(tOnOffBlocks); %number of on and off phases
            
            tBlockIdx = 1; %the index of phase of that trial (phase = tBlock)           
   
            previousSpike = start; %questionable
            spikeIdx = 1;
            lastY = 0;
            
            %spikeYvalues{iChann, iTrial}(1) = 0;
            
            while tBlockIdx < numTblock
               
                    currentSpike = consideredSpikes{iChann, iTrial}(spikeIdx);
                    if previousSpike > tOnOffTimes(tBlockIdx + 1)
                        tBlockIdx = tBlockIdx + 1;
                    end
                    % This fist if statement accounts for when the prevoius
                    % and current spikes are in the same on or off phase.
                    % In this case, the boolean caseBool will be true
                    if tBlockIdx <= numTblock && currentSpike <= stop && previousSpike >= tOnOffTimes(tBlockIdx) && currentSpike <= tOnOffTimes(tBlockIdx +1)
                        
                        [spikeYvalues{iChann, iTrial}(spikeIdx), lastY, onOrOff] = findSpikeYval2(firingRatesOn(iChann), firingRatesOff(iChann), previousSpike, currentSpike,tOnOffTimes, tOnOffBlocks, tBlockIdx, lastY, 1);
                       
                    % This else statement accounts for when the previous 
                    % spike and current spike are in different on or off 
                    % phases. In this case, the boolean caseBool is false.             
                    else
                        if tBlockIdx <= numTblock && currentSpike <= stop && previousSpike >= tOnOffTimes(tBlockIdx) && previousSpike <= tOnOffTimes(tBlockIdx + 1) && currentSpike >= tOnOffTimes(tBlockIdx +1)
                            [spikeYvalues{iChann, iTrial}(spikeIdx), lastY, onOrOff] = findSpikeYval2(firingRatesOn(iChann), firingRatesOff(iChann), previousSpike, currentSpike,tOnOffTimes, tOnOffBlocks, tBlockIdx, lastY, 0);
                            
                        end
                    end
                previousSpike = currentSpike;
                spikeIdx = spikeIdx + 1;    
               
            end
        end
    end
end
