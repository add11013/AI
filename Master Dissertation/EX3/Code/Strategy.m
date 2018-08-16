for iiii=1:10
OpenOriginalData=[];
OriginalData=[];
load(['Data_EX2' int2str(iiii)])

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

for i=1:length(PredictClose)
    RO=RealOpen(i,1);
    RC=RealClose(i,1);
    PC=PredictClose(i,1);
    if RO>PC
        sale=sale+1;
        if RO>RC
            EarnMoney=EarnMoney+(RO-RC);
        else
            LostMoney=LostMoney+(RC-RO);
        end
    end
    if RO<PC
        buy=buy+1;
        if RO<RC
            EarnMoney=EarnMoney+(RC-RO);
        else
            LostMoney=LostMoney+(RO-RC);
        end
    end
end
Profit=EarnMoney-LostMoney;
AAAA(iiii,1)=Profit;
AAAA(iiii,2)=Sale;
AAAA(iiii,3)=buy;
end