function [Profit,SellAndBuy]=PaperStrategy(Actual,Forecast)
output_ProfitTable=[];
tempProfitTable=[];

for iii=1:300
    Alpha=iii*0.001-0.001;
    buy=0;
    sell=0;
    SaleMoney=0;
    BuyMoney=0;
    for i=1:length(Forecast)-1
        %sell
        if abs(Forecast(i)-Actual(i))/Actual(i)<=Alpha && (Forecast(i+1)-Actual(i))>0
            sell=sell+1;
            SaleMoney=SaleMoney+Actual(i+1)-Actual(i);
        end
        %buy
        if abs(Forecast(i)-Actual(i))/Actual(i)<=Alpha && (Forecast(i+1)-Actual(i))<0
            buy=buy+1;
            BuyMoney=BuyMoney+Actual(i)-Actual(i+1);
        end
    end
    Profit(iii,1)=SaleMoney+BuyMoney;
    SellAndBuy(iii,1)=sell;
    SellAndBuy(iii,2)=buy;
end
end