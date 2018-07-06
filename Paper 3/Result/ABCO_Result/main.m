clear
clc
close all

NumberOfTrial=10;
LoadName='ABCOResult_EX3_trial';

%% calculate the profit and the operation of the training data
ProfitTable=0;
OperateTable=0;
for Trial=1:NumberOfTrial
        load([LoadName int2str(Trial)])
        % every taget have to calculate the profit

        TMP=[];
        j1=1;
        % if number of target is one, get the first output real part
        if NumberOfTarget==1
            for N=1:NumberOfOUTPUT
                RealValue=OriginalData(32:(31+NumberOfTrainPoint),j1);
                ModelOutput=real(ABCOgBest.yHead(:,N))+OriginalData(31:(31+NumberOfTrainPoint-1),j1);
            end
        else
            for N=1:NumberOfOUTPUT
                RealValue(:,j1)=OriginalData(32:(31+NumberOfTrainPoint),j1);
                RealValue(:,j1+1)=OriginalData(32:(31+NumberOfTrainPoint),j1);
                ModelOutput(:,j1)=real(ABCOgBest.yHead(:,N))+OriginalData(31:(31+NumberOfTrainPoint-1),j1);
                ModelOutput(:,j1+1)=imag(ABCOgBest.yHead(:,N))+OriginalData(31:(31+NumberOfTrainPoint-1),j1);
                j1=j1+2;
            end
        end
    for Target=1:NumberOfTarget
        Actual=RealValue(:,Target);
        Forecast=ModelOutput(:,Target);
        Interval=30;
        Differential=OriginalData((32-Interval):NumberOfTrainPoint+31);
        for i=1:length(Differential)-Interval
            tmp=[];
            for ii=1:Interval
                tmp(ii)=Differential(i+Interval-ii+1)-Differential(i+Interval-ii);
            end
            IntervalMean(i)=sum(tmp)/Interval;
            IntervalStd(i)=std(tmp);
        end
        
        % In function CalculateProfit_1, the final input parameter decides
        % whether find the best alpha, if 0 means yes
        [TrainProfit(Target).t(Trial).value,TrainOperationTable(Target).t(Trial).value]=CalculateProfit_2(Actual,Forecast,IntervalMean,IntervalStd,0);
    end
end


%% caluculate the sum of the all targets profit
for Trial=1:NumberOfTrial
        TrainAllProfit(Trial).value=0;
        for Target=1:NumberOfTarget
            TrainAllProfit(Trial).value=TrainAllProfit(Trial).value+TrainProfit(Target).t(Trial).value;
        end
end
save('TrainProfitTable','TrainAllProfit','TrainProfit','TrainOperationTable','NumberOfTrial','LoadName','NumberOfTrial','Interval');
clear
load('TrainProfitTable');

%% calculate the mean of the all trials for training data
for Trial=1:NumberOfTrial
    for i=1:100
        TrainAllProfitTable(i,Trial)=TrainAllProfit(Trial).value(i,1);
    end
end
        for i=1:100
            Mean(i,1)=mean(TrainAllProfitTable(i,:));
        end
%% find the best alpha        
[Maximum MaxIndex]=max(Mean);
MaxIndex=(MaxIndex+1)*0.001;
Maximum
['The best threshold parameter in the training data is: ' num2str(MaxIndex)]
BestAlpha=MaxIndex;

%% plot the profits in different alpha
% x=linspace(1,100,100);
% hold on 
% plot(x,Mean)

%% use best alpha to calculate the profits of the test data
ProfitTable=0;
OperateTable=0;
for Trial=1:NumberOfTrial
        load([LoadName int2str(Trial)])
        % every taget have to calculate the profit

        TMP=[];
        j1=1;
        % if number of target is one, get the first output real part
        if NumberOfTarget==1
            for N=1:NumberOfOUTPUT
                RealValue(:,j1)=real(OriginalData((32+NumberOfTrainPoint):(31+NumberOfAllPoint),N));
                ModelOutput(:,j1)=real(testyHead(:,N))+OriginalData((31+NumberOfTrainPoint):(31+NumberOfAllPoint-1),j1);
            end
        else
            for N=1:NumberOfOUTPUT
                RealValue(:,j1)=OriginalData((32+NumberOfTrainPoint):(31+NumberOfAllPoint),j1);
                RealValue(:,j1+1)=OriginalData((32+NumberOfTrainPoint):(31+NumberOfAllPoint),j1);
                ModelOutput(:,j1)=real(testyHead(:,N))+OriginalData((31+NumberOfTrainPoint):(31+NumberOfAllPoint-1),j1);
                ModelOutput(:,j1+1)=imag(testyHead(:,N))+OriginalData((31+NumberOfTrainPoint):(31+NumberOfAllPoint-1),j1);
                j1=j1+2;
            end
        end
    for Target=1:NumberOfTarget
        Actual=RealValue(:,Target);
        Forecast=ModelOutput(:,Target);
        Differential=OriginalData((32+NumberOfTrainPoint-Interval):NumberOfAllPoint+31);
        for i=1:length(Differential)-Interval
            tmp=[];
            for ii=1:Interval
                tmp(ii)=Differential(i+Interval-ii+1)-Differential(i+Interval-ii);
            end
            IntervalMean(i)=sum(tmp)/Interval;
            IntervalStd(i)=std(tmp);
        end
        % in function CalculateProfit_1, the final input parameter decides
        % whether find the best alpha, if 0 means yes
        [Profit(Target).t(Trial).value,OperationTable(Target).t(Trial).value]=CalculateProfit_2(Actual,Forecast,IntervalMean,IntervalStd,0);
    end
end

%% caluculate the sum of the all targets profit
for Trial=1:NumberOfTrial
        AllProfit(Trial).value=0;
        for Target=1:NumberOfTarget
            AllProfit(Trial).value=AllProfit(Trial).value+Profit(Target).t(Trial).value;
            TargetProft(Target).value(:,Trial)=Profit(Target).t(Trial).value;
        end
end
save('ProfitTable','AllProfit','Profit','OperationTable','NumberOfTrial','LoadName','NumberOfTrial','BestAlpha','TargetProft');
clear
load('ProfitTable');

%% calculate the mean of the all trials for training data
for Trial=1:NumberOfTrial
    for i=1:100
        AllProfitTable(i,Trial)=AllProfit(Trial).value(i,1);
    end
end
for i=1:100
            Mean(i)=mean(AllProfitTable(i,:));
            Std(i)=std(AllProfitTable(i,:));
end
BestAlpha=roundn(BestAlpha*1000,-2);
['The mean of the all trials profits which by the best threshold parameter: ' int2str(Mean(BestAlpha))]
['The standard deviation of the all trials profits which by the best threshold parameter: ' int2str(Std(BestAlpha))]

%MaxValue_IN_BestAlpha is the max profit in the best training threshold
%parameter
%Max_MODEl is the model which the max profit is
[MaxValue_IN_BestAlpha Max_MODEL]=max(AllProfitTable(BestAlpha,:));

%% find the best alpha
%[Maximum MaxIndex]=max(Mean);
[SortingValue SortingIndex]=sort(AllProfitTable(:,Max_MODEL));
SortingValue=flipud(roundn(SortingValue,-2));
SortingIndex=flipud(roundn(SortingIndex,-5)*0.001);