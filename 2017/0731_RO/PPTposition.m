j1=1;
for i=1:4
    j2=1;
    for jj=1:3
        A(j2:j2+3,i)=PSO.gBestPosition(j1:j1+3);
        j1=j1+4;
        j2=j2+4;
    end
end