% This is for running paper experimentations
for trial=2:10
    
    trial
    ExperimentationName='E6';
    NumberOfTrainPoint=115;
    ParameterOfPreSC=0.3;
    ParameterOfConsSC=0.3;
    MODEL_PSO(trial, ExperimentationName, NumberOfTrainPoint, ParameterOfPreSC, ParameterOfConsSC);
end