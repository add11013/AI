figure(5)
xx=linspace(1,NumberOfAllPoint,NumberOfAllPoint);
temp=real(ACOgBest.Error(:,1))+real(ACOgBest.Error(:,2))+imag(ACOgBest.Error(:,1))+imag(ACOgBest.Error(:,2));
temp=temp/4;
temp1=real(testError(:,1))+real(testError(:,2))+imag(testError(:,1))+imag(testError(:,2));
temp1=temp1/4;
temp=[temp; temp1];
plot(xx,temp)
line([NumberOfTrainPoint NumberOfTrainPoint],[min(real(ACOgBest.Error(:,1)))*1.5 max(real(ACOgBest.Error(:,1)))*1.5]);
