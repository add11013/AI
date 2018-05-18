for iiii=1:10
    j=iiii;
    chr = int2str(j);
    load(['PSOResult_EX1_trial' chr])
    aaaTrain(iiii,1)=PSOgBest.Distance;
    temp0=0;
       for N=1:NumberOfOUTPUT
            temp=sqrt(sum(testError(:,N).*conj(testError(:,N)))/(NumberOfTrainPoint));
            temp0=temp0+temp;
       end
       aaaTrain(iiii,2)=temp0;
end