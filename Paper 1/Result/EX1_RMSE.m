for iiii=1:10
    j=iiii;
    chr = int2str(j);
    load(['Result_E1_trial' chr]);
    rmse(iiii,:)=[RMSENASO2 RMSENASC2 RMSESPO2 RMSESPC2];
end