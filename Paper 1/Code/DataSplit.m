% This function is for data splitting
%parameter "data" is Original data
function output=DataSplit(ExperimentationName)
    data=xlsread(['data_' ExperimentationName '.csv']);
    data=data(1:size(data,1),:);
    for i=1:size(data,2)
        A=data(:,i);
        B(:,i)=A(~isnan(A));
    end
    output=B;
end