function ProfitTable=PaperStrategy(LoadName,Target,Alpha,NumberOfTest)
for iiii=1:NumberOfTest
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
    
    
    
    actual=OriginalData(33:length(OriginalData),:);
    buy=0;
    sale=0;
    SaleMoney=0;
    BuyMoney=0;
    Forecast=PredictClose(:,Target);
    for i=1:length(Forecast)-1
        %sell
        if abs(Forecast(i)-actual(i))/actual(i)<=Alpha && (Forecast(i+1)-actual(i))>0
            sale=sale+1;
            SaleMoney=SaleMoney+actual(i+1)-actual(i);
        end
        
        if abs(Forecast(i)-actual(i))/actual(i)<=Alpha && (Forecast(i+1)-actual(i))<0
            buy=buy+1;
            BuyMoney=BuyMoney+actual(i)-actual(i+1);
        end
    end
    Profit=SaleMoney+BuyMoney;
    ProfitTable(iiii,1)=Profit;
    ProfitTable(iiii,2)=sale;
    ProfitTable(iiii,3)=buy;
end
end