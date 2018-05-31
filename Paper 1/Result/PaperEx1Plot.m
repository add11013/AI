figure(5)
xx=linspace(1,NumberOfAllPoint,NumberOfAllPoint);
temp=real(PSOgBest.Error(:,1))+real(PSOgBest.Error(:,2))+imag(PSOgBest.Error(:,1))+imag(PSOgBest.Error(:,2));
temp=temp/4;
temp1=real(testError(:,1))+real(testError(:,2))+imag(testError(:,1))+imag(testError(:,2));
temp1=temp1/4;
temp=[temp; temp1];
plot(xx,temp)
line([NumberOfTrainPoint NumberOfTrainPoint],[min(real(PSOgBest.Error(:,1)))*1.5 max(real(PSOgBest.Error(:,1)))*1.5]);

x=linspace(1,length(OriginalData)-31,length(OriginalData)-31);

        %畫出每個目標
        for i=1:NumberOfTarget
            figure(i)
            hold on
            plot(x,OriginalData(32:length(OriginalData),i));
        end
        
        x=linspace(1,NumberOfAllPoint,NumberOfAllPoint);
        %目標為第31天的漲跌，所以 第31天價格+第31天漲跌(預測出的)=第32天價格(預測出的)
        for N=1:NumberOfOUTPUT
            FinalyHead(:,N)=[PSOgBest.yHead(:,N) ;testyHead(:,N)];
        end
        k=1;

        figure(1)
        output=OriginalData(31:length(OriginalData)-1,1)+real(FinalyHead(:,1));
        plot(x,output);
        line([NumberOfTrainPoint NumberOfTrainPoint],[min(OriginalData(:,1))*0.5 max(OriginalData(:,1))*1.3]);
        legend('Target (NASDAQ Opening)','Forecast (NASDAQ Opening)');
        xlabel('Trading date index');
        ylabel('Stock price (NASDAQ Opening)');
        
        figure(2)
        output=OriginalData(31:length(OriginalData)-1,2)+imag(FinalyHead(:,1));
        plot(x,output);
        line([NumberOfTrainPoint NumberOfTrainPoint],[min(OriginalData(:,2))*0.8 max(OriginalData(:,2))*1.3]);
        legend('Target (NASDAQ Closing)','Forecast (NASDAQ Closing)');
        xlabel('Trading date index');
        ylabel('Stock price (NASDAQ Closing)');
        
        figure(3)
        output=OriginalData(31:length(OriginalData)-1,3)+real(FinalyHead(:,2));
        plot(x,output);
        line([NumberOfTrainPoint NumberOfTrainPoint],[min(OriginalData(:,3))*0.8 max(OriginalData(:,3))*1.3]);
        legend('Target (S&P500 Opening)','Forecast (S&P500 Opening)');
        xlabel('Trading date index');
        ylabel('Stock price (S&P500 Opening)');
        
        figure(4)
        output=OriginalData(31:length(OriginalData)-1,4)+imag(FinalyHead(:,2));
        plot(x,output);
        line([NumberOfTrainPoint NumberOfTrainPoint],[min(OriginalData(:,4))*0.8 max(OriginalData(:,4))*1.3]);
        legend('Target (S&P500 Closing)','Forecast (S&P500 Closing)');
        xlabel('Trading date index');
        ylabel('Stock price (S&P500 Closing)');
        
        figure(7)
        plot(1:PSO.iterations,PSO.plotRMSE);
        legend('Learning curve');
        xlabel('Iteration');
        ylabel('RMSE');