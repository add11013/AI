clear
clc
 
Target=1;
NumberOfTest=10;
Alpha=0.07;
LoadName='PSOResult_Ex2_trial';
LoadName1='1PSOResult_Ex2_trial';

ProfitTable=[PaperStrategy(LoadName,Target,Alpha,NumberOfTest) PaperStrategy(LoadName1,Target,Alpha,NumberOfTest)];

Mean(1)=mean(ProfitTable(:,1));
Mean(2)=mean(ProfitTable(:,4));