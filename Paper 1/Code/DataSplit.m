% This function is for data splitting
%parameter "data" is Original data
function output=DataSplit(ExperimentationName)
    data=xlsread(['data_' ExperimentationName '.csv']);
    data=data(1:size(data,1),:);
    for i=1:size(data,2)
        A=data(:,i);
        B(:,i)=A(~isnan(A));
    end
    %target 6
    if strcmp(ExperimentationName,'E6')
        for i=1:size(data,2)
            output(:,i)=log(B(2:length(B),i)./B(1:length(B)-1,i));
        end
    else
        output=B;
    end
end