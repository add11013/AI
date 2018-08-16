
for j=1:10
    load(['Result_E5_trial' int2str(j)])
    aRMSETable(j)=PSOgBest.Distance;
end
