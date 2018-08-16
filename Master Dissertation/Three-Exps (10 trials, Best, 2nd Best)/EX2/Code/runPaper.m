% This is for running paper experimentations
NumberOfTrial=10;
% for trial=1:NumberOfTrial
%     trial
%     ExperimentationName='EX1';
%     NumberOfTrainPoint=181;
%     ParameterOfPreSC=0.2;
%     MODEL_PSO1(trial, ExperimentationName, NumberOfTrainPoint, ParameterOfPreSC);
% end
% 
for trial=1:NumberOfTrial
    trial
    ExperimentationName='EX2';
    NumberOfTrainPoint=204;
    ParameterOfPreSC=0.3;
    MODEL_PSO(trial, ExperimentationName, NumberOfTrainPoint, ParameterOfPreSC);
end
% 
for trial=1:NumberOfTrial
    trial
    ExperimentationName='EX3';
    NumberOfTrainPoint=181;
    ParameterOfPreSC=0.2;   
    MODEL_PSO(trial, ExperimentationName, NumberOfTrainPoint, ParameterOfPreSC);
end
% 
% for trial=1:NumberOfTrial
%     trial
%     ExperimentationName='EX4';
%     NumberOfTrainPoint=800;
%     ParameterOfPreSC=0.1;
%     MODEL_PSO(trial, ExperimentationName, NumberOfTrainPoint, ParameterOfPreSC);
% end
% 
% for trial=1:NumberOfTrial
%     trial
%     ExperimentationName='EX5';
%     NumberOfTrainPoint=400;
%     ParameterOfPreSC=0.2;
%     MODEL_PSO(trial, ExperimentationName, NumberOfTrainPoint, ParameterOfPreSC);
% end

% for trial=1:NumberOfTrial
%     trial
%     ExperimentationName='EX6';
%     NumberOfTrainPoint=115;
%     ParameterOfPreSC=0.3;
%     MODEL_PSO(trial, ExperimentationName, NumberOfTrainPoint, ParameterOfPreSC);
% end


