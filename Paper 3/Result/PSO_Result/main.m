clear
clc
 
Target=1;
NumberOfTest=10;

LoadName='PSOResult_Ex3_trial';
LoadName1='1PSOResult_Ex3_trial';

ProfitTable=[RunStrategy(LoadName,Target,NumberOfTest) RunStrategy(LoadName1,Target,NumberOfTest)];

Mean(1)=mean(ProfitTable(:,1));
Mean(2)=mean(ProfitTable(:,4));