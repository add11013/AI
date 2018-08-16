for ii=1:4
    for mm=1:121
        for jj=1:121
            if isnan(IIM(ii).value(mm,jj))
                 IIM(ii).value(mm,jj)=0;
            end
        end
    end
end
