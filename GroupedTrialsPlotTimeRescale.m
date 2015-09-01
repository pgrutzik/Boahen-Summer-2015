function [groupedTrialsZ] = GroupedTrialsPlotTimeRescale(~)

clear all;
set(0, 'defaulttextinterpreter', 'latex');

load('~/Documents/BoahenLabProject15/poolAllData_numState_2.mat');
load('~/Documents/BoahenLabProject15/hash.mat');

[z, tao] = timeRescaleISI(popData, popParam);

allIdx = strcmp(popParam.recording, '2011_12_17');
usedConditions = find(allIdx);
start = usedConditions(1);
stop = usedConditions(end);

groupedTrialsZ = cell(length(popData.states),1); % MAKE SURE LENGHT MATCHES NUMBER OF TRIALS TESTED
sortedZs = cell(length(popData.states),1);
b = cell(length(popData.states),1);
confidence = cell(length(popData.states),1);

%for n =  start : stop 
    figure;
 n = 257;
    %Loop through all of the channels 
    sortedZs{n} = cell(16,1);
    b{n} = cell(16,1);
    confidence{n} = cell(16,1);
    for c = 1:16
        
        % Assume there are at least 100 spikes in each trial. 
        numSpikes = 100;
        
        % The first column of groupedTrials is the b values and the second
        % column is the sorted z values
        groupedTrialsZ{n}{c} = zeros(length(popParam.behaviorParam{n}.trialNumbers) * numSpikes,1);
        sortedZs{n}{c} = zeros(length(popParam.behaviorParam{n}.trialNumbers) * numSpikes,1);
        
        % You have to keep track of the last index of the nonZero values of
        % groupedTrialsZsort and groupedTrialsBvals so that you know where
        % to start imputting values of the next trial
        
        %Loop through all fo the trials and sort all of the trials create a
        %very large vector of all of the z values. Then loop through all of
        %those z values and add the b values to the second column of the
        %groupedTrials matrix.
        
        lastIndex = 1;
        for y = 1:length(popParam.behaviorParam{n}.trialNumbers) % for each trial
        
            l_z = length(z{n}{y}{c});
            groupedTrialsZ{n}{c}(lastIndex:lastIndex + (l_z-1)) = z{n}{y}{c};  %Maybe I should just add all of the z's and then sort later

            lastIndex = lastIndex + l_z;
        end 
        groupedTrialsZ{n}{c} = groupedTrialsZ{n}{c}(groupedTrialsZ{n}{c} ~= 0);
        sortedZs{n}{c} = sort(groupedTrialsZ{n}{c});
        b{n}{c}  = zeros(length(sortedZs{n}{c}),1);
        confidence{n}{c} = zeros(length(sortedZs{n}{c}),1);
        %Loop through aand create b vector
        for x = 1:length(sortedZs{n}{c})
            b{n}{c}(x) = (x - 0.5)/length(sortedZs{n}{c});
            confidence{n}{c}(x,1) = b{n}{c}(x) - (1.36 * length(sortedZs{n}{c}))^ (-0.5);
            confidence{n}{c}(x,2) = b{n}{c}(x) + (1.36 * length(sortedZs{n}{c}))^ (-0.5);
        end
%         figure;
%         hist(b{n}{c});
%         xlabel('b values');
%         ylabel('Density');
%         title('Probability Density of b values');
%         print('Probability Density of b values', '-dpdf');
%         
%         figure(1);
%         hist(sortedZs{n}{c},100);
%         xlabel('z values');
%         ylabel('Density');
%         title('Probability Density of z values');
%         print('Probability Density of z values', '-dpdf');
%         
        hold on;
        hist(tao{n}{c},100);
        xlabel('tao values');
        ylabel('Density');
        title('Probability Density of tao values');
        print('Probability Density of tao values', '-dpdf');
        
%         figure;
%         cdfplot(sortedZs{n}{c});
%         xlabel('z values');
%         ylabel('Cumulative Density');
%         title('Cumulative Density of z values');
%         print('Cumulative Density of z values', '-dpdf');
%         
%         figure;
%         cdfplot(b{n}{c});
%         xlabel('b values');
%         ylabel('Cumulative Density');
%         title('Cumulative Density of b values');
%         print('Cumulative Density of b values', '-dpdf');
        
%         figure;
%         hold on;
%         plot(b{n}{c}, sortedZs{n}{c}, 'b')
%         
%         lengthConfidence = length(confidence{n}{c});
%         spikes = 0:1/lengthConfidence: 1-(1/lengthConfidence);
%         axis([0,1,0,1]);
%         xlabel('b values');
%         ylabel('z values');
%         title('Kolmogorov Smirnov (KS) Plot');
%         plot(spikes, confidence{n}{c}(:,1), 'r')
%         plot(spikes, confidence{n}{c}(:,2), 'r')
%         print('KS Plot', '-dpdf');
%         hold off;
    end
%end
        
end