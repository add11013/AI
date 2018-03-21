%% testing part
        %NumberOfTarget指的是實數型態的目標個數
        NumberOfTarget=size(OriginalData,2);
        NumberOfAllPoint=size(DataMatrix,1);
        NumberOfTestPoint=NumberOfAllPoint-NumberOfTrainPoint;
        % get  Testing data
        %get h1~hM
        for M=1:length(FeatureIndex)
            testh(M).value=DataMatrix(NumberOfTrainPoint+1:NumberOfAllPoint,FeatureIndex(M));
        end
        
        %get test target
        k=61;
        for N=1:NumberOfOUTPUT
            realPartOfTest=DataMatrix(NumberOfTrainPoint+1:NumberOfAllPoint,k);
            imagPartOfTest=DataMatrix(NumberOfTrainPoint+1:NumberOfAllPoint,k+1);
            k=k+2;
            testy(N).value=realPartOfTest+imagPartOfTest*j;
        end
            
        %把gBest放進termSet、Lambda1Set、Lambda2Set的集合
          %將粒子的位子，儲存成termSet，以便後續建造前艦部的fuzzy set
            j1=1;
            for M=1:NumberOfINPUT
                for number=1:length(h(M).center)
                    termSet.INPUT(M).fuzzyset(number).value=[PSOgBest.Position(j1) PSOgBest.Position(j1+1)];
                    Lambda1Set.INPUT(M).fuzzyset(number)=PSOgBest.Position(j1+2);
                    Lambda2Set.INPUT(M).fuzzyset(number)=PSOgBest.Position(j1+3);
                    j1=j1+4;
                end
            end

                
        %每個OUTPUT都要算Beta
        for N=1:NumberOfOUTPUT
            %算每一條規則的啟動強度
            for jj=1:NumberOfTestPoint
                for rule=1:length(formationMatrix)
                    membership=1;
                    for M=1:NumberOfINPUT
                        r(M)=gaussmf(testh(M).value(jj),termSet.INPUT(M).fuzzyset(formationMatrix(rule,M)).value,1);
                        theata1Ofh=gaussmf(testh(M).value(jj),termSet.INPUT(M).fuzzyset(formationMatrix(rule,M)).value,3)*Lambda1Set.INPUT(M).fuzzyset(formationMatrix(rule,M)); %dr/dx
                        theata2Ofh=gaussmf(testh(M).value(jj),termSet.INPUT(M).fuzzyset(formationMatrix(rule,M)).value,6)*Lambda2Set.INPUT(M).fuzzyset(formationMatrix(rule,M)); %d2r/dx^2
                        temp=r(M)*exp(j*(theata1Ofh+theata2Ofh));
                        membership=membership*temp;
                    end
                    testBeta(N).value(rule,jj)=membership;
                end
            end
        end
        
        
        for N=1:NumberOfOUTPUT
            %把每個beta裡的規則分開儲存
            for K=1:NumberOfPremise
                B(N).k(K).value=testBeta(N).value(K,:);
            end
        end


      %算Q個P
        for Q=1:NumberOfCons
            for jj=1:NumberOfTestPoint
                %CaculateB
                BOfLambdaTMP0=[];
                for N=1:NumberOfOUTPUT
                    TMP0=[];
                    for K=1:NumberOfPremise
                        TMP=[B(N).k(K).value(jj)];
                        TMP0=[TMP0 TMP];
                    end
                    BOfLambdaTMP=[TMP0];
                    BOfLambdaTMP0=[BOfLambdaTMP0; BOfLambdaTMP];
                end
                BOfLambda=BOfLambdaTMP0;
                
                %CaculateL
                LOfLambdaTMP0=[];
                for K=1:NumberOfPremise
                    TMP0=[];
                    for N=1:NumberOfOUTPUT
                        TMP=[lambda(N).k(K).q(Q).value(jj)];
                        TMP0=[TMP0 TMP];
                    end
                    LOfLambdaTMP=[TMP0];
                    LOfLambdaTMP0=[LOfLambdaTMP0;LOfLambdaTMP];
                end
                LOfLambda=LOfLambdaTMP0;
                
                %CaculateH
                HOfLambdaTMP0=[];
                for N=1:NumberOfOUTPUT
                    TMP0=[1];
                    %TMP0=[1 h1 h2 ... hM]
                    for M=1:NumberOfINPUT
                        TMP=testh(M).value(jj);
                        TMP0=[TMP0 TMP];
                    end
                    HOfLambdaTMP=[TMP0];
                    HOfLambdaTMP0=[HOfLambdaTMP0;HOfLambdaTMP];
                end
                HOfLambda=HOfLambdaTMP0;
                
                testP.q(Q).value(jj,:)=[BOfLambda*LOfLambda*HOfLambda];
            end
        end
        
        %
        TMP0=[];
        for Q=1:NumberOfCons
            TMP=testP.q(Q).value;
            TMP0=[TMP0 TMP];
        end
        testA=TMP0;

      for jj=1:NumberOfTestPoint
        for N=1:NumberOfOUTPUT
            yyy=testA(jj,:)*PSOgBest.Theata;  %testyHead
            testyHead(jj,N)=yyy(N,1);
           %caculate error
            testError(jj,N)=testy(N).value(jj,1)-testyHead(jj,N);
        end
      end