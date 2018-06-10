clear
clc
 
Target=1;
NumberOfTrial=10;
LoadName='ABCOResult_EX3_trial';
LoadName1='1ABCOResult_Ex3_trial';
ProfitTable=PaperStrategy(LoadName,Target,NumberOfTrial);
ProfitTable1=PaperStrategy(LoadName1,Target,NumberOfTrial);
for i=1:300
    Mean(i)=mean(ProfitTable(i).value(:,1));
    Mean1(i)=mean(ProfitTable1(i).value(:,1));
end

% ProfitTable=[PaperStrategy(LoadName,Target,NumberOfTrial) PaperStrategy(LoadName1,Target,NumberOfTrial)];

% Mean(StrategyAlpha+1)=mean(ProfitTable(:,1));
% Mean(StrategyAlpha+1)=mean(ProfitTable(:,4));