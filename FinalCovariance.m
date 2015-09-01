function [meanCovar,covH, covP, meanCorrCoef, corrH, corrP, onOffPairs, corrRecH,corrRecP] = FinalCovariance(popParam, popData_states, popData_timeBins, TonVToffVBoth)
set(0, 'defaulttextinterpreter', 'latex');
%start = 2;
%popData should be an array of cells where each element is an array of
%cells. Each of these cells should have an array:   { {[   ]} {[     ]} {[
%]} }
% 
%TonVToffVBoth is a string input either 
%           'onOff' - makes the program compare only the effect of on 
%           phase times on off phase times
%           'offOn' - makes the program compare only the effect of off 
%           phase times on on phase times
%           'both' - makes the program compare the effect of off 
%           phase times on on phase times and on phase times on off phase
%           times

%set(0, 'defaulttextinterpreter', 'latex');
binSize = 0.01;
numBins = cell(length(popParam.iCond),1);
avgNumBins = zeros(length(popParam.iCond), 1);
covar = zeros(length(popParam.iCond), 1);
corrCOEF = zeros(length(popParam.iCond), 1);
corrCoef = zeros(length(popParam.iCond), 1);
onOffPairs = cell(length(popParam.iCond),1); %this num for now, think it grows as you go
%avgToffCondition = zeros(608, 1);
%avgTonCondition = zeros(608, 1);
for conditionIdx = 1:length(popParam.iCond)  %numel(attendIdx)
   
%for n = 1: length(popData_states) %looping through conditions (recording x orientation x attention, 608)
  % n = 257;
    %TonToff = zeros(length(popData_states{conditionIdx}), 1);
    numBins{conditionIdx} = NaN(length(popData_states{conditionIdx}),1);
    totalToff = 0;
    totalTon = 0;
    countTon = 0;
    countToff = 0;
    a = 1;
    %calculate largest number of bins you could have. 
    max_length = sum(cellfun(@length, popData_states{conditionIdx}));
    onOffPairs{conditionIdx} = zeros(max_length,2);
    for s = 1:length(popData_states{conditionIdx}) % loop through trials of each condition
        %new_a = a;
        numBins{conditionIdx}(s) = length(popData_states{conditionIdx}{s});
        i = 1;
        %Reset Toff and Ton at the start of every trial 
        Toff = 0;
        Ton = 0;

        while i < numBins{conditionIdx}(s) %loop through the number of time bins
            
            % find Ton times and add them together
            if(popData_states{conditionIdx}{s}(i) == 2)
                onStartIndex = i;
                while (i < numBins{conditionIdx}(s) && popData_states{conditionIdx}{s}(i) == 2)
                    i = i + 1;
                end
                onStopIndex = i; % -1
                if i == numBins{conditionIdx}(s) + 1
                    Ton = popData_timeBins{conditionIdx}{s}(end) - popData_timeBins{conditionIdx}{s}(onStartIndex) + binSize;
                else
                    Ton = popData_timeBins{conditionIdx}{s}(onStopIndex) - popData_timeBins{conditionIdx}{s}(onStartIndex); % one Ton time. should we count them up seperately or can you just add them to a vector of some sort and call mean(vec)?
                end
                totalTon = totalTon + Ton;
                countTon = countTon +1;
            end
            if ~strcmp(TonVToffVBoth,'onOff')
                if (Toff ~= 0)
                    %TonToff(a) = Ton * Toff;
                    onOffPairs{conditionIdx}(a,1) = Ton;
                    onOffPairs{conditionIdx}(a,2) = Toff;
                    a = a + 1;
                end
            end
            if (i < numBins{conditionIdx}(s) && popData_states{conditionIdx}{s}(i) == 1)
                % find Toff times and add them together
                offStartIndex = i;
                while (i < numBins{conditionIdx}(s) && popData_states{conditionIdx}{s}(i) == 1)
                    i = i + 1;
                end
                offStopIndex = i; % -1;
                if i == numBins{conditionIdx}(s) +1
                    Toff = popData_timeBins{conditionIdx}{s}(end) - popData_timeBins{conditionIdx}{s}(offStartIndex) + binSize;
                else
                    Toff = popData_timeBins{conditionIdx}{s}(offStopIndex) - popData_timeBins{conditionIdx}{s}(offStartIndex); % one Toff time
                end
                totalToff = totalToff + Toff;
                countToff = countToff +1;
            end
            %TonVToffVBoth SKIP 
            if ~strcmp(TonVToffVBoth,'offOn') 
                if (Ton ~= 0)
                    %TonToff(a) = Ton * Toff;
                    onOffPairs{conditionIdx}(a,1) = Ton;
                    onOffPairs{conditionIdx}(a,2) = Toff;
                    a = a + 1;
                end
            end
        end %loop through time Bins
        % at the end: take average of Toff times
        
        
    end % loop through trials
    %trims all of the unset elements in onOffPairs
    onOffPairs{conditionIdx}(a:end,:) = [];
    
    avgTon2 = mean(onOffPairs{conditionIdx}(:,1));
    avgToff2 = mean(onOffPairs{conditionIdx}(:,2));
    TonToff2 = zeros(length(onOffPairs{conditionIdx}(:,1)),1);
    for m = 1: size(onOffPairs{conditionIdx}(:,1))
        TonToff2(m) = onOffPairs{conditionIdx}(m,1) * onOffPairs{conditionIdx}(m,2);
    end
    meanTonToff2 = mean(TonToff2);
    covar(conditionIdx) = meanTonToff2 - (avgTon2 * avgToff2);
    divisor = sqrt(nanvar(onOffPairs{conditionIdx}(:,1))* nanvar(onOffPairs{conditionIdx}(:,2)));
    corrCOEF(conditionIdx) = covar(conditionIdx)/divisor;
    one = onOffPairs{conditionIdx}(:,1);
    two = onOffPairs{conditionIdx}(:,2);
    %corrCoef(conditionIdx) = corrcoef(one, two); 
    avgNumBins(conditionIdx) = nanmean(numBins{conditionIdx});
    
end % loop through conditions

%loop through conditions and separate by recording
%grab all of the unique recordings
recordings = unique(popParam.recording);
conditionByRec = NaN(length(recordings),1);

for r = 1: length(recordings)
    Rrec = strcmp(popParam.recording(:),recordings(r));
    idxOfR = find(Rrec);
    numOfRec = sum(Rrec);
    conditionByRec(r) = nanmean(corrCOEF(idxOfR(1):idxOfR(end)));        
end

    
    [corrRecH,corrRecP] = ttest(conditionByRec);
    figure;
    hold on;
    m = 1: length(recordings);
    z = 1:length(avgNumBins);
    plot(m, conditionByRec, 'k+')
    %plot(z, avgNumBins, 'r')
    xlabel('Unique Recording');
    ylabel('Correlation Coefficient');
    title('Correlation Coefficient of Length of On and Off States Grouped By Recording');
    %legend('Average Number of Time Bins in Trials', 'Correlation Coefficient');
    print('Off_On Correlation Coefficient by Recording', '-dpdf');
    hold off;

    meanCovar = nanmean(covar);
    [covH,covP] = ttest(covar); % H: 1 for false, 0 for true
    meanCorrCoef = nanmean(corrCOEF);
    [corrH,corrP] = ttest(corrCOEF); % H: 1 for false, 0 for true
    
    endCOEF = nanmean(corrCOEF(400:end))
    
    figure;
    hold on;
    y = 1:numel(corrCOEF);
    subplot(2,1,1);
    plot(y,corrCOEF,'+')
    xlabel('Condition Number');
    ylabel('Correlation Coefficient');
    title('Correlation Coefficient of Length of On and Off States');   
    subplot(2,1,2);
    plot(z, avgNumBins, 'r')
    xlabel('Condition Number');
    ylabel('Average Number of Time Bins in Trials');
    print('OffOn Correlation Coefficient all Conditions', '-dpdf');
    hold off;
    
    x = 1:numel(covar);
    figure;
    hold on;
    subplot(2,1,1);
    plot(x, covar, 'o')  
    xlabel('Condition Number');
    ylabel('Covariance (s2)');
    title('Covariance of Length of On and Off States ');
    subplot(2,1,2);
    plot(z, avgNumBins, 'r') 
    xlabel('Condition Number');
    ylabel('Average Number of Time Bins in Trials');
    print('OffOn Covariance all Conditions', '-dpdf');
    hold off;
end