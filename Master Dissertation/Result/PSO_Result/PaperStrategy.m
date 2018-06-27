function [Profit,operation]=PaperStrategy(Actual,Forecast,IntervalMean,IntervalStd,Alpha)
%if Alpha ==0 means that we have to find the best alpha, so we change the
%Alpha in every for loop
        if Alpha==0
            flag=true;
        else
            flag=false;
        end

output_ProfitTable=[];
tempProfitTable=[];

for i=1:100
    if flag==true
        Alpha=i*0.001;
    end
    buy=0;
    sell=0;
    SaleMoney=0;
    BuyMoney=0;
    for t=1:length(Forecast)-1
        %buy
        if abs(Forecast(t)-Actual(t))/Actual(t)<=Alpha && (Forecast(t+1)-Actual(t))>0
            operation(t,i)=1;
            BuyMoney=BuyMoney+Actual(t+1)-Actual(t);
        %sell
        elseif abs(Forecast(t)-Actual(t))/Actual(t)<=Alpha && (Forecast(t+1)-Actual(t))<0
            operation(t,i)=-1;
            SaleMoney=SaleMoney+Actual(t)-Actual(t+1);
        else
            operation(t,i)=0;
        end
    end
    Profit(i,1)=SaleMoney+BuyMoney;
end
end