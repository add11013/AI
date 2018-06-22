function [Profit,operation]=PaperStrategy(Actual,Forecast,TestRMSE)
output_ProfitTable=[];
tempProfitTable=[];

for iii=1:100
    Alpha=iii*0.001-0.001;
    Alpha1=Alpha-0.03;
    buy=0;
    sell=0;
    SaleMoney=0;
    BuyMoney=0;
    for i=1:length(Forecast)-1
        %buy
        if abs(Forecast(i)-Actual(i))/Actual(i)<=Alpha && abs(Forecast(i)-Actual(i))/Actual(i)>Alpha1 && (Forecast(i+1)-Actual(i))>0
            operation(i,iii)=1;
            BuyMoney=BuyMoney+Actual(i+1)-Actual(i);
        %sell
        elseif abs(Forecast(i)-Actual(i))/Actual(i)<=Alpha && abs(Forecast(i)-Actual(i))/Actual(i)>Alpha1 && (Forecast(i+1)-Actual(i))<0
            operation(i,iii)=-1;
            SaleMoney=SaleMoney+Actual(i)-Actual(i+1);
        else
            operation(i,iii)=0;
        end
    end
    Profit(iii,1)=SaleMoney+BuyMoney;
%     SellAndBuy(iii,1)=sell;
%     SellAndBuy(iii,2)=buy;
end
end