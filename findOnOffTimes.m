function [stopStartTimes] = findOnOffTimes(popData)
load('~/Documents/BoahenLabProject15/poolAllData_numState_2.mat');
load('~/Documents/BoahenLabProject15/hash.mat');
stopStartTimes = cell(length(popData.states),1);
%[stopStartTimes{ = cell(length(popData.states),1); [X{1:3, 1}] = deal(zeros(3))
for n = 1: length(popData.states) %looping through conditions (recording x orientation x attention, 608)
%     a = 2;
stopStartTimes{n} = cell(length(popData.states{n}),1);
    for s = 1:length(popData.states{n}) % loop through trials of each condition
        
        numBins = length(popData.states{n}{s}); %number of time bins in trial
        stopStartTimes{n}{s} = zeros(numBins,2);
        i = 1;
        startIndex = 1;
        while i < numBins % i = the specific 10 ms time bin
            % find Ton times and add them together
            if(popData.states{n}{s}(i) == 1)
               
                stopStartTimes{n}{s}(startIndex,1) = popData.timeBins{n}{s}(i); %SO INSTEAD OF MAKING THIS INDEX AT i MAKE IT AT ANOTHER SPECIFIED INDEX
                stopStartTimes{n}{s}(startIndex,2) = 1;
                startIndex = startIndex + 1;
                while (i <= numBins && popData.states{n}{s}(i) == 1)
                    i = i + 1;
                end
                if i-1 == numBins 
                    % just changed the (onStartIndex) to (i)
                    stopStartTimes{n}{s}(startIndex,1) = popData.timeBins{n}{s}(i-1)+0.01; %i think this actually just rewrites the same number in the same place
                    stopStartTimes{n}{s}(startIndex,2) = 0; 
                    i = i + 1; %DONT THINK YOU NEED THIS
                    
                end
            end

            if (i < numBins && popData.states{n}{s}(i) == 2)
                % find Toff times and add them together
                offStartIndex = i; %WHY DO YOU HAVE + 0.01 ONLY ON THE OFF STATES? JUST DELETED + 0.01 FROM LINE BELOW
                stopStartTimes{n}{s}(startIndex,1) = popData.timeBins{n}{s}(i); %i think this actually just rewrites the same number in the same place
                stopStartTimes{n}{s}(startIndex,2) = 2;
                startIndex = startIndex + 1;
                while (i <= numBins && popData.states{n}{s}(i) == 2)
                    i = i + 1;
                end
                if i-1 == numBins 
                    offStopIndex = i; % JUST CHANGED (OFFSTARTINDEX) TO (I-1)
                    stopStartTimes{n}{s}(startIndex,1) = popData.timeBins{n}{s}(i-1) + 0.01; %i think this actually just rewrites the same number in the same place
                    stopStartTimes{n}{s}(startIndex,2) = 0;
                    i = i + 1; %DONT THINK YOU NEED THIS
                  
                end
            end
        end
        stopStartTimes{n}{s}( all(~stopStartTimes{n}{s},2), : ) = [];
    end % loop through trials

end

end