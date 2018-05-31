
for j=1:10
    j
    name=['Result_E1_trial' int2str(j) ];
    load(name);
%     swarm=rmfield(swarm,'RLSE');
    clear MAPEAPPLE1
    clear MAPEAPPLE2
    clear MAPEDELL1
    clear MAPEDELL2
    clear MAPEIBM1
    clear MAPEIBM1
    clear MAPEMicroSoft1
    clear MAPEMicroSoft2

    save(name);
end