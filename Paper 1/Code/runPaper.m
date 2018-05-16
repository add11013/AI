% This is for running paper experimentations
for trial=1:1
    
    trial
    ExperimentationName='E3';
    NumberOfTrainPoint=433;
    ParameterOfPreSC=0.15;
    ParameterOfConsSC=0.07;
    MODEL_PSO(trial, ExperimentationName, NumberOfTrainPoint, ParameterOfPreSC, ParameterOfConsSC);
end