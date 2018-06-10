clear
clc
 
Target=2;
NumberOfTest=10;
Alpha=0.07;
LoadName='PSOResult_Ex3_trial';
LoadName1='1PSOResult_Ex3_trial';

ProfitTable=[PaperStrategy(LoadName,Target,Alpha,NumberOfTest) PaperStrategy(LoadName1,Target,Alpha,NumberOfTest)];

Mean(1)=mean(ProfitTable(:,1));
Mean(2)=mean(ProfitTable(:,4));