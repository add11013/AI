function ProfitTable=RunStrategy(LoadName,Target,NumberOfTest);
for iiii=1:NumberOfTest
OpenOriginalData=[];
OriginalData=[];
load([LoadName int2str(iiii)])

for N=1:NumberOfOUTPUT
            FinalyHead(:,N)=[PSOgBest.yHead(:,N) ;testyHead(:,N)];
end
TMP=[];
j1=1;
for N=1:NumberOfOUTPUT
    TMP=OriginalData(31:length(OriginalData)-1,j1)+real(FinalyHead(:,N));
    PredictClose(:,j1)=TMP(1:NumberOfAllPoint-1);
    TMP=OriginalData(31:length(OriginalData)-1,j1+1)+imag(FinalyHead(:,N));
    PredictClose(:,j1+1)=TMP(1:NumberOfAllPoint-1);
    j1=j1+2;
end


    RealClose=OriginalData(33:length(OriginalData),:);
    RealOpen=OpenOriginalData(33:length(OriginalData),:);
    
    buy=0;
    sale=0;
    EarnMoney=0;
    LostMoney=0;
    Lost=0;
    Earn=0;

for i=1:length(PredictClose)
    RO=RealOpen(i,Target);
    RC=RealClose(i,Target);
    PC=PredictClose(i,Target);
    if RO>PC
        sale=sale+1;
        if RO>RC
            EarnMoney=EarnMoney+(RO-RC);
            Earn=Earn+1;
        else
            LostMoney=LostMoney+(RC-RO);
            Lost=Lost+1;
        end
    end
    if RO<PC
        buy=buy+1;
        if RO<RC
            EarnMoney=EarnMoney+(RC-RO);
            Earn=Earn+1;
        else
            LostMoney=LostMoney+(RO-RC);
            Lost=Lost+1;
        end
    end
end
Profit=EarnMoney-LostMoney;
ProfitTable(iiii,1)=Profit;
ProfitTable(iiii,2)=sale;
ProfitTable(iiii,3)=buy;
end
