function [output_ProfitTable]=PaperStrategy(LoadName,Target,NumberOfTrial)
output_ProfitTable=[];
tempProfitTable=[];
for iiii=1:NumberOfTrial
    OriginalData=[];
    load([LoadName int2str(iiii)])
    for iii=1:300
        Alpha=iii*0.001-0.001
        
        for N=1:NumberOfOUTPUT
            FinalyHead(:,N)=[ABCOgBest.yHead(:,N) ;testyHead(:,N)];
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
        
        
        
        actual=OriginalData(33:length(OriginalData),Target);
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
        ProfitTable(1)=Profit;
        ProfitTable(2)=sale;
        ProfitTable(3)=buy;
        tempProfitTable=[tempProfitTable;ProfitTable];
        output_ProfitTable(iii).value(iiii,:)=ProfitTable;
    end
end

end