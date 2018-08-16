% % This is for running paper experimentations
for trial=1:10
    trial
    ExperimentationName='EX1';
    NumberOfTrainPoint=181;
    ParameterOfPreSC=0.2;
    MODEL_ABCO(trial, ExperimentationName, NumberOfTrainPoint, ParameterOfPreSC);
end

for trial=1:10
    trial
    ExperimentationName='EX2';
    NumberOfTrainPoint=204;
    ParameterOfPreSC=0.3;
    MODEL_ABCO(trial, ExperimentationName, NumberOfTrainPoint, ParameterOfPreSC);
end

for trial=1:10
    trial
    ExperimentationName='EX3';
    NumberOfTrainPoint=400;
    ParameterOfPreSC=0.15;
    MODEL_ABCO(trial, ExperimentationName, NumberOfTrainPoint, ParameterOfPreSC);
end
