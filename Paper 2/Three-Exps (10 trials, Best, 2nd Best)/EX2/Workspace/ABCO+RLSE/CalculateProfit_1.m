function [ Profit,OperationTable ] = CalculateProfit_1( Actual,Forecast,IntervalMean,IntervalStd,Alpha)
%if Alpha ==0 means that we have to find the best alpha, so we change the
%Alpha in every for loop
        if Alpha==0
            flag=true;
        else
            flag=false;
        end
        %calculate profit in differnt alpha
        for i=1:100
            if flag==true
                Alpha=i*0.001;
            end
            SellMoney=0;
            BuyMoney=0;
            for t=1:length(Forecast)-1
                Threshold=(IntervalMean(t)+IntervalStd(t)*0);
                % when forecast value is bigger than today value => buy
                if abs(Forecast(t)-Actual(t))/Actual(t)<=Alpha && (Forecast(t+1)-Actual(t))>0
                    if (Forecast(t+1)-Actual(t))>Threshold && Threshold>0
                        OperationTable(t,i)=1;
                        BuyMoney=BuyMoney+Actual(t+1)-Actual(t);
                    else
                        OperationTable(t,i)=0;
                    end
                % when forecast value is smaller than today value => sell
                elseif abs(Forecast(t)-Actual(t))/Actual(t)<=Alpha && (Forecast(t+1)-Actual(t))<0
                    if (Actual(t)-Forecast(t+1))>abs(Threshold) && Threshold<0
                        OperationTable(t,i)=-1;
                        SellMoney=SellMoney+Actual(t)-Actual(t+1);
                    else
                        OperationTable(t,i)=0;
                    end
                else
                    OperationTable(t,i)=0;
                end
            end
            Profit(i,1)=SellMoney+BuyMoney;
        end
end

