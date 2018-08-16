function Select=Roulette(P)

m = length(P);
Select = zeros(1,1);
r = rand(1,1);

    sumP = 0;
    j = ceil(m*rand);    %?生1~m之?的?机整?
    while sumP < r
        sumP = sumP + P(mod(j-1,m)+1);
        j = j+1;
    end
    Select = mod(j-2,m)+1;
end
