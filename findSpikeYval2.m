function [y, lastY, onOrOff] = findSpikeYval2(firingRateOn, firingRateOff, previousSpike, currentSpike, tOnOffTimes, tOnOffBlocks, tBlockIdx, lastY, caseBool)
   % do I have way too many parameters to be acceptable? 
   % if onOrOff = 0, ISI is during off position. if onOrOff = 1, ISI is
   % during on position
    ISI = currentSpike - previousSpike; % maybe previous spike is being tampered with?
    if caseBool == 1       
        if tOnOffBlocks(tBlockIdx) == 1
            onOrOff = 0;
            y = (firingRateOff * ISI) + lastY;
        else 
            onOrOff = 1;
            y = (firingRateOn * ISI) + lastY;          
        end
    else
        firstChunk = tOnOffTimes(tBlockIdx + 1) - previousSpike;
        secondChunk = currentSpike - tOnOffTimes(tBlockIdx +1);
               
        if tOnOffBlocks(tBlockIdx) == 1
            onOrOff = 0;
            firstY = (firingRateOff * firstChunk) + lastY;
            y = (firingRateOn * secondChunk) + firstY;
        else
            onOrOff = 1;
            firstY = (firingRateOn * firstChunk) + lastY;
            y = (firingRateOff * secondChunk) + firstY;
        end
    end
    lastY = y;   
end