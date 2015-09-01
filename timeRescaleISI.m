function [z, tao] = timeRescaleISI(popData, popParam)

[stopStartTimes] = findOnOffTimes(popData);

[spikeYvalues, ISI] = rescaledSpikeValues(stopStartTimes); %this is cell of all conditions, cell of all channels, cell of all trails, y values

% alln = strcmp(popParam.recording, '2011_12_17');
% consideredRecordings = find(alln);
% start = startIndex(1);
% stop = startIndex(end);
% finalSSTimes = cell(length(startIndex),1);

allIdx = strcmp(popParam.recording, '2011_12_17');
usedConditions = find(allIdx);
tao = cell(length(popData.states),1);
z = cell(length(popData.states),1);
start = usedConditions(1);
stop = usedConditions(end);

%for n =  start : stop 
n = 257;
%start : stop  % Loop through all condtitions that are in the recording you want (GO BACK AND FIGURE THIS OUT AT THE END)
    
    tao{n} = cell(length(popParam.behaviorParam{n}.trialNumbers),1);
    z{n} = cell(length(popParam.behaviorParam{n}.trialNumbers),1);
    
    for y = 1:length(popParam.behaviorParam{n}.trialNumbers) % for each trial
        %figure;
        tao{n}{y} = cell(16,1);
        z{n}{y} = cell(16,1);
        
        for c = 1:16  % for each channel a(a~=0)
            %notZeroIndexY = find(spikeYvalues{n}{y}{c});
            channelspikeY = spikeYvalues{n}{y}{c}(spikeYvalues{n}{y}{c} ~= 0); %MAYBE GO BACK AND PUT THIS INTO RESCALEDSPIKEVALUES.M
            length_channelspikeY = length(channelspikeY);
            tao{n}{y}{c} = zeros(length_channelspikeY-1,1);
            z{n}{y}{c} = zeros(length_channelspikeY-1,1);
            if ~isempty(spikeYvalues{n}{y}{c})
                %tao{n}{y}{c}(1) = spikeYvalues{n}{y}{c}(1);
                for q = 1:length_channelspikeY-1
                    tao{n}{y}{c}(q) = channelspikeY(q+1)- channelspikeY(q); % These are the differences between the sequentail rescaled spike times
                    %           this means tao{n}{y}{c}(1)
                    %             z{n}{y}{c}(q) = 1 - exp(-tao{n}{y}{c}(q));
                    %             b{n}{y}{c}(q) = (q - 0.5)/(length_channelsikeY-1);
                end
                length_tao = length(tao{n}{y}{c});
                nonZero_tao = tao{n}{y}{c}(tao{n}{y}{c} ~= 0); % currently between 100 and 99 because the index
                length_nonZero_tao = length(nonZero_tao);
                z{n}{y}{c} = NaN(length_nonZero_tao,1);  %WHY DO I INITIALIZE Z AND B TWICE???

                for p = 1:length_nonZero_tao % this different approach causes slight but not significant changes in b
                    exponent = tao{n}{y}{c}(p);
                    z{n}{y}{c}(p) = 1 - exp(-(exponent));

                end
            end
            
%             hold on;
%             hist(tao{n}{y}{c},100);
%             xlabel('Tao values');
%             ylabel('Density');
%             title('Probability Density of Tao values');
%             print('Probability Density of Tao values', '-dpdf');
            % make an array of the average distance of points on the graph away from
            % the 45 degree line. Then plot this line. You will be able to see
            % how and by how much these points increase.
            
            
            %now plot all of the trials of each channel on the same plot. So
            %there will be 16 figures. I am going to have to do this in a
            %separate function. I have to export all of the zSort and the b
        end
    end
%end
end

