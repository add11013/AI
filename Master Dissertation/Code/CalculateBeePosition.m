function  NewPosition = CalculateBeePosition(RandStart,RandEnd,max,min)

    %-3~5  =>  8*rand()-3

    randtemp=(RandStart+RandEnd)*rand()+RandStart;
    NewPosition=min+randtemp*(max-min);

end

