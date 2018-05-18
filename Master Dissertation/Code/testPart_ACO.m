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
        k=30*NumberOfTarget+1;
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
                    termSet.INPUT(M).fuzzyset(number).value=[ACOgBest.Position(j1) ACOgBest.Position(j1+1)];
                    Lambda1Set.INPUT(M).fuzzyset(number)=ACOgBest.Position(j1+2);
                    Lambda2Set.INPUT(M).fuzzyset(number)=ACOgBest.Position(j1+3);
                    j1=j1+4;
                end
            end
            
         
         %每個OUTPUT都要算Beta
            %算每一條規則的啟動強度

                for rule=1:length(FormationMatrix)
                    testmembership1=1;
                    testmembership2=1;
                    testmembership3=1;
                    for M=1:NumberOfINPUT
                        r=gaussmf(testh(M).value,termSet.INPUT(M).fuzzyset(FormationMatrix(rule,M)).value,1);
                        theata1Ofh=gaussmf(testh(M).value,termSet.INPUT(M).fuzzyset(FormationMatrix(rule,M)).value,3)*Lambda1Set.INPUT(M).fuzzyset(FormationMatrix(rule,M)); %dr/dx
                        theata2Ofh=gaussmf(testh(M).value,termSet.INPUT(M).fuzzyset(FormationMatrix(rule,M)).value,6)*Lambda2Set.INPUT(M).fuzzyset(FormationMatrix(rule,M)); %d2r/dx^2
                        temp=r.*exp(j.*(theata1Ofh+theata2Ofh));
                        testmembership1=testmembership1.*temp;
                        temp2=r.*cos(theata2Ofh).*cos(theata1Ofh)+r.*cos(theata2Ofh).*sin(theata1Ofh).*j;
                        testmembership2=testmembership2.*temp2;
%                         temp3=r*cos(theata2Ofh)*sin(theata1Ofh)+r*sin(theata2Ofh)*j;
%                         membership3=membership3*temp3;
                    end
                    testBeta(1).value(rule,:)=testmembership1;
                    testBeta(2).value(rule,:)=testmembership2;
%                     testBeta(3).value(rule,jj)=membership3;
                end

            %Normalization
            for N=1:NumberOfOUTPUT
                for jj=1:NumberOfTestPoint
                    for rule=1:length(FormationMatrix)
                        temp1=real(testBeta(N).value(rule,jj))/sum(real(testBeta(N).value(:,jj)));
                        temp2=imag(testBeta(N).value(rule,jj))/sum(imag(testBeta(N).value(:,jj)));
                        testNBeta(N).value(rule,jj)=temp1+temp2*j;
                    end
                end
            end
            %replace beta with normalizationed value 
                testBeta=testNBeta;
        
        
        for N=1:NumberOfOUTPUT
            %把每個beta裡的規則分開儲存
            for K=1:NumberOfPremise
                testB(N).k(K).value=testBeta(N).value(K,:);
            end
        end
        
        temp=[];        
        for jj=1:NumberOfTestPoint      
            temp2=[];
            for N=1:NumberOfOUTPUT
                temp1=[];
                for K=1:NumberOfPremise    
                    testlambda=testB(N).k(K).value(jj);
                    temp=testlambda;
                    for M=1:NumberOfINPUT
                        temp=[temp testlambda*h(M).value(jj)];
                    end
                    temp1=[temp1 temp];
                end
                temp2=[temp2; temp1];                
            end
            testA(jj).value=temp2;
        end

        
        %new_yHead(output)
        for jj=1:NumberOfTestPoint
             testyHead(jj,:)=testA(jj).value*ACOgBest.theta;  %y
        end
        %caculate error    
       for N=1:NumberOfOUTPUT
            testError(:,N)=testy(N).value-testyHead(:,N);
       end

            
            