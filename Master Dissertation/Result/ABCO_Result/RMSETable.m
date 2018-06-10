clear
for iiii=1:10
    j=iiii;
    chr = int2str(j);
    load(['ABCOResult_EX2_trial' chr])
    aaaTrain(iiii,1)=ABCOgBest.rmse;
    temp0=0;
       for N=1:NumberOfOUTPUT
            temp=sqrt(testError(:,N)'*(testError(:,N))/(NumberOfTestPoint));
            temp0=temp0+temp;
       end
       aaaTrain(iiii,2)=temp0;
end
