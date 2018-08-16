function weight = CalculateWeight(q, k, l)
%CALCULATEWEIGHT Summary of this function goes here
%   Detailed explanation goes here
gg=gaussmf(l,[1 q*k],1);
weight=1./(q.*k.*sqrt(2*pi))*gg;
end

