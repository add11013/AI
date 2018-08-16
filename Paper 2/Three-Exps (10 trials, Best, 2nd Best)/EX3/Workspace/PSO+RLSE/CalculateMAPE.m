for iiii=5:5
    j=iiii;
    chr = int2str(j);
    load(['Result_EX3_trial' chr]);
    temp0=[];
    temp1=[];
    temp0=OriginalData(32:length(OriginalData),1);
    temp1=OriginalData(31:length(OriginalData)-1,1)+real(FinalyHead(:,1));
    temp0_1=temp0(1:NumberOfTrainPoint);
    temp1_1=temp1(1:NumberOfTrainPoint);
    temp0_2=temp0(NumberOfTrainPoint+1:NumberOfAllPoint);
    temp1_2=temp1(NumberOfTrainPoint+1:NumberOfAllPoint);
    
    MAPE1_1=sum(abs(temp0_1-temp1_1)./temp0_1)/NumberOfTrainPoint*100;
    MAPE1_2=sum(abs(temp0_2-temp1_2)./temp0_2)/NumberOfTestPoint*100;
    
    AAAAA(iiii,1)=MAPE1_2;
    
    temp0=[];
    temp1=[];
    temp0=OriginalData(32:length(OriginalData),2);
    temp1=OriginalData(31:length(OriginalData)-1,2)+imag(FinalyHead(:,1));
    temp0_1=temp0(1:NumberOfTrainPoint);
    temp1_1=temp1(1:NumberOfTrainPoint);
    temp0_2=temp0(NumberOfTrainPoint+1:NumberOfAllPoint);
    temp1_2=temp1(NumberOfTrainPoint+1:NumberOfAllPoint);
    
    MAPE2_1=sum(abs(temp0_1-temp1_1)./temp0_1)/NumberOfTrainPoint*100;
    MAPE2_2=sum(abs(temp0_2-temp1_2)./temp0_2)/NumberOfTestPoint*100;
    
    AAAAA(iiii,2)=MAPE2_2;
    
    

end