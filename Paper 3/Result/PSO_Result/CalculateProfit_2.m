    function [ Profit,OperationTable ] = CalculateProfit_2( Actual,Forecast,IntervalMean,IntervalStd,Alpha)
    %if Alpha ==0 means that we have to find the best alpha, so we change the
    %Alpha in every for loop
    if Alpha==0
        flag=true;
    else
        flag=false;
    end
    WindowInterval=5;
    Windows=length(Actual)/WindowInterval;
    for i=1:100
        t=1;
        HoldStock=0;
        if flag==true
            Alpha=i*0.001;
        end
        Money=0;

        for ii=1:Windows-1
            for iii=1:WindowInterval
                Threshold=IntervalMean(t)+IntervalStd(t)*Alpha;
                %buy
                if abs(Forecast(t)-Actual(t))/Actual(t)<=Alpha && (Forecast(t+1)-Actual(t))>0
                    if (Forecast(t+1)-Actual(t))>Threshold && Threshold>0
                        OperationTable(t,i)=1;
                        HoldStock=HoldStock+1;
                        Money=Money-Actual(t);
                    else
                        OperationTable(t,i)=0;
                    end
                    %sell
                elseif HoldStock>0 && abs(Forecast(t)-Actual(t))/Actual(t)<=Alpha && (Actual(t)-Forecast(t+1))<0
                    if (Actual(t)-Forecast(t+1))>Threshold && Threshold<0
                        OperationTable(t,i)=-1;
                        HoldStock=HoldStock-1;
                        Money=Money+Actual(t);
                    else
                        OperationTable(t,i)=0;
                    end

                    %do nothing
                else
                    OperationTable(t,i)=0;
                end
                t=t+1;
            end
            if HoldStock~=0
                Money=Money+HoldStock*Actual(t);
                HoldStock=0;
            end
        end
        % calculate the final window
        for iii=1:mod(length(Actual),WindowInterval)
             Threshold=IntervalMean(t)+IntervalStd(t)*0;
            %buy
            if abs(Forecast(t)-Actual(t))/Actual(t)<=Alpha && (Forecast(t+1)-Actual(t))>0
                if (Forecast(t+1)-Actual(t))>Threshold && Threshold>0
                    OperationTable(t,i)=1;
                    HoldStock=HoldStock+1;
                    Money=Money-Actual(t);
                else
                    OperationTable(t,i)=0;
                end
                %sell
            elseif HoldStock>0 && abs(Forecast(t)-Actual(t))/Actual(t)<=Alpha && (Actual(t)-Forecast(t+1))<0
                if (Actual(t)-Forecast(t+1))>abs(Threshold) && Threshold<0
                    OperationTable(t,i)=-1;
                    HoldStock=HoldStock-1;
                    Money=Money+Actual(t);
                else
                    OperationTable(t,i)=0;
                end

                %do nothing
            else
                OperationTable(t,i)=0;
            end
            t=t+1;
        end
        if HoldStock~=0
            Money=Money+HoldStock*Actual(t);
            HoldStock=0;
        end

        Profit(i,1)=Money;
    end
    end

