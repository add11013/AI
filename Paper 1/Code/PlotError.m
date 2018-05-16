xx=linspace(1,NumberOfAllPoint,NumberOfAllPoint);
temp=[PSOgBest.Error; testError]
plot(xx,temp)
line([NumberOfTrainPoint NumberOfTrainPoint],[min(PSOgBest.Error(:,1))*0.8 max(PSOgBest.Error(:,1))*1.3]);
