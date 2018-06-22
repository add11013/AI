clear
clc
close all

NumberOfProfitTarget=2;
NumberOfTrial=10;
LoadName='PSOResult_EX2_trial';

%% 100 iterations
ProfitTable=0;
OperateTable=0;
for Target=1:NumberOfProfitTarget
    for Trial=1:NumberOfTrial
        load([LoadName int2str(Trial)])
        for N=1:NumberOfOUTPUT
            FinalyHead(:,N)=[PSOgBest.yHead(:,N) ;testyHead(:,N)];
        end
        TMP=[];
        j1=1;
        if NumberOfProfitTarget==1
            for N=1:NumberOfOUTPUT
                TMP=OriginalData(31:length(OriginalData)-1,j1)+real(FinalyHead(:,N));
                PredictClose(:,j1)=TMP(1:NumberOfAllPoint-1);
            end
        else
            for N=1:NumberOfOUTPUT
                TMP=OriginalData(31:length(OriginalData)-1,j1)+real(FinalyHead(:,N));
                PredictClose(:,j1)=TMP(1:NumberOfAllPoint-1);
                TMP=OriginalData(31:length(OriginalData)-1,j1+1)+imag(FinalyHead(:,N));
                PredictClose(:,j1+1)=TMP(1:NumberOfAllPoint-1);
                j1=j1+2;
            end
        end

        
        
        
        Actual=OriginalData(33:length(OriginalData),Target);        
        Forecast=PredictClose(:,Target);
        TestRMSE=PSOgBest.Distance;
        [Profit,OperationTable]=PaperStrategy(Actual,Forecast,TestRMSE);
        
        Profit_tmp(:,Trial)=Profit;
    end
    ProfitTable=ProfitTable+Profit_tmp;
    
end
save('ProfitTable','ProfitTable','OperationTable','NumberOfProfitTarget','NumberOfTrial','LoadName');
clear
load('ProfitTable','NumberOfProfitTarget','NumberOfTrial','LoadName');
%% 1 iteration
LoadName=['1' LoadName];
one_iteration_ProfitTable=0;
one_iteration_OperateTable=0;
for Target=1:NumberOfProfitTarget
for Trial=1:NumberOfTrial
        load([LoadName int2str(Trial)])
        for N=1:NumberOfOUTPUT
            FinalyHead(:,N)=[PSOgBest.yHead(:,N) ;testyHead(:,N)];
        end
        TMP=[];
        j1=1;
        if NumberOfProfitTarget==1
            for N=1:NumberOfOUTPUT
                TMP=OriginalData(31:length(OriginalData)-1,j1)+real(FinalyHead(:,N));
                PredictClose(:,j1)=TMP(1:NumberOfAllPoint-1);
            end
        else
            for N=1:NumberOfOUTPUT
                TMP=OriginalData(31:length(OriginalData)-1,j1)+real(FinalyHead(:,N));
                PredictClose(:,j1)=TMP(1:NumberOfAllPoint-1);
                TMP=OriginalData(31:length(OriginalData)-1,j1+1)+imag(FinalyHead(:,N));
                PredictClose(:,j1+1)=TMP(1:NumberOfAllPoint-1);
                j1=j1+2;
            end
        end
        
        
        
        Actual=OriginalData(33:length(OriginalData),Target);
        Forecast=PredictClose(:,Target);
        TestRMSE=PSOgBest.Distance;
        [Profit,one_iteration_OperationTable]=PaperStrategy(Actual,Forecast,TestRMSE);
        Profit_tmp(:,Trial)=Profit;
end
        one_iteration_ProfitTable=one_iteration_ProfitTable+Profit_tmp;
end
save('one_iteration_ProfitTable','one_iteration_ProfitTable','one_iteration_OperationTable');
clear

%% mean
load('ProfitTable')
load('one_iteration_ProfitTable')
for i=1:100
        Mean(i,1)=mean(ProfitTable(i,:));
        Std(i,1)=std(ProfitTable(i,:));
        Mean(i,2)=mean(one_iteration_ProfitTable(i,:));
        Std(i,2)=std(one_iteration_ProfitTable(i,:));
end
max(Mean(:,1))
max(Mean(:,2))

x=linspace(1,100,100);
hold on 
plot(x,Mean(:,1))
plot(x,Mean(:,2),'--')