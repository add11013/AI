
for j=1:10
    load(['Result_E4_trial' int2str(j)])
    aRMSETable(j)=PSOgBest.Distance;
end
