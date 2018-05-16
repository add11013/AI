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
        k=3*NumberOfTarget+1;
        for N=1:NumberOfOUTPUT
            realPartOfTest=DataMatrix(NumberOfTrainPoint+1:NumberOfAllPoint,k);
            imagPartOfTest=DataMatrix(NumberOfTrainPoint+1:NumberOfAllPoint,k+1);
            k=k+2;
            testy(N).value=realPartOfTest+imagPartOfTest*j;
        end
            
          %將粒子的位置，儲存成termSet，以便後續建造前艦部的fuzzy set
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
            %算每一條規則的啟動強度
            for jj=1:NumberOfTestPoint
                for rule=1:length(FormationMatrix)
                    testmembership1=1;
                    testmembership2=1;
                    testmembership3=1;
                    for M=1:NumberOfINPUT
                        r=gaussmf(testh(M).value(jj),termSet.INPUT(M).fuzzyset(FormationMatrix(rule,M)).value,1);
                        theata1Ofh=gaussmf(testh(M).value(jj),termSet.INPUT(M).fuzzyset(FormationMatrix(rule,M)).value,3)*Lambda1Set.INPUT(M).fuzzyset(FormationMatrix(rule,M)); %dr/dx
                        theata2Ofh=gaussmf(testh(M).value(jj),termSet.INPUT(M).fuzzyset(FormationMatrix(rule,M)).value,6)*Lambda2Set.INPUT(M).fuzzyset(FormationMatrix(rule,M)); %d2r/dx^2
                        temp=r*exp(j*(theata1Ofh+theata2Ofh));
                        testmembership1=testmembership1*temp;
                        temp2=r*cos(theata2Ofh)*cos(theata1Ofh)+r*cos(theata2Ofh)*sin(theata1Ofh)*j;
                        testmembership2=testmembership2*temp2;
%                         temp3=r*cos(theata2Ofh)*sin(theata1Ofh)+r*sin(theata2Ofh)*j;
%                         membership3=membership3*temp3;
                    end
                    testBeta(1).value(rule,jj)=testmembership1;
                    testBeta(2).value(rule,jj)=testmembership2;
%                     testBeta(3).value(rule,jj)=membership3;
                end
            end
            %Normalization
            for N=1:2
                for jj=1:NumberOfTestPoint
                    for rule=1:length(FormationMatrix)
                        temp1=real(testBeta(N).value(rule,jj))/sum(real(testBeta(N).value(:,jj)));
                        temp2=imag(testBeta(N).value(rule,jj))/sum(imag(testBeta(N).value(:,jj)));
                        testNBeta(N).value(rule,jj)=temp1+temp2*j;
                    end
                end
            end
            %replace beta with normalizationed value 
            for N=1:2
                testBeta(N).value=testNBeta(N).value;
            end
        
        
        for N=1:NumberOfOUTPUT
            %把每個beta裡的規則分開儲存
            for K=1:NumberOfPremise
                testB(N).k(K).value=testBeta(N).value(K,:);
            end
        end
        
        %轉換箭靶中心
        for N=1:NumberOfOUTPUT
            for Q=1:NumberOfCons
                    for K=1:NumberOfPremise
                        r1=exp(-(testB(N).k(K).value-aocC(N).q(Q)).*conj(testB(N).k(K).value-aocC(N).q(Q))./(2.*aocS(N).q(Q).*conj(aocS(N).q(Q))));
                        w1=r1.*real(-(testB(N).k(K).value-aocC(N).q(Q))/aocS(N).q(Q).*conj(aocS(N).q(Q)));
                        testlambda(N).k(K).q(Q).value=r1.*exp(j*w1);
                    end
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
                        TMP=[testB(N).k(K).value(jj)];
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
                        TMP=[testlambda(N).k(K).q(Q).value(jj)];
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
                        TMP=h(M).value(jj);
                        TMP0=[TMP0 TMP];
                    end
                    HOfLambdaTMP=[TMP0];
                    HOfLambdaTMP0=[HOfLambdaTMP0;HOfLambdaTMP];
                end
                HOfLambda=HOfLambdaTMP0;
                
                testP.q(Q).point(jj).value=BOfLambda*LOfLambda*HOfLambda;
            end
        end
        
        
        for N=1:NumberOfOUTPUT
            TMP0=[];
            TMP=[];
            for Q=1:NumberOfCons
                for jj=1:NumberOfTestPoint
                    TMP(jj,:)=testP.q(Q).point(jj).value(N,:);
                end
                TMP0=[TMP0 TMP];
            end
            testRLSE(N).A=TMP0;
        end


          
      %new_yHead(output)
        for N=1:NumberOfOUTPUT
            testyHead(:,N)=testRLSE(N).A*PSOgBest.RLSE.Theata;  %y
           %caculate error
            testError(:,N)=testy(N).value-testyHead(:,N);
        end
