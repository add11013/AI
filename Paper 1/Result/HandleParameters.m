
for i=1:10
    i
    name=['Result_E3_trial' int2str(i) ];
    load(name);
    swarm=rmfield(swarm,'RLSE');
    swarm=rmfield(swarm,'A');
    clear lambda;
    clear testlambda;
    save(name);
end