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
        if NumberOfTarget==1
            testy(1).value=DataMatrix(NumberOfTrainPoint+1:NumberOfAllPoint,k);
        elseif mod(NumberOfTarget,2)==0
            for N=1:NumberOfOUTPUT
                realPartOfTest=DataMatrix(NumberOfTrainPoint+1:NumberOfAllPoint,k);
                imagPartOfTest=DataMatrix(NumberOfTrainPoint+1:NumberOfAllPoint,k+1);
                k=k+2;
                testy(N).value=realPartOfTest+imagPartOfTest*j;
            end
        else
            for N=1:NumberOfOUTPUT
                realPartOfTest=DataMatrix(NumberOfTrainPoint+1:NumberOfAllPoint,k);
                imagPartOfTest=DataMatrix(NumberOfTrainPoint+1:NumberOfAllPoint,k+1);
                k=k+2;
                testy(N).value=realPartOfTest+imagPartOfTest*j;
            end
            testy(N+1).value=DataMatrix(NumberOfTrainPoint+1:NumberOfAllPoint,k);
        end
            
          %將粒子的位置，儲存成termSet，以便後續建造前艦部的fuzzy set
            j1=1;
            for M=1:NumberOfINPUT
                for number=1:length(h(M).center)
                    termSet.INPUT(M).fuzzyset(number).value=[PSOgBest.Position(j1) PSOgBest.Position(j1+1)];
                    Lambda1Set.INPUT(M).fuzzyset(number)=PSOgBest.Position(j1+2);
                    j1=j1+3;
                end
            end
            
         
         %每個OUTPUT都要算Beta
            %算每一條規則的啟動強度

                for rule=1:length(FormationMatrix)
                    testmembership1=1;
                    testmembership2=1;
                    for M=1:NumberOfINPUT
                        r=gaussmf(testh(M).value,termSet.INPUT(M).fuzzyset(FormationMatrix(rule,M)).value,1);
                        theata1Ofh=gaussmf(testh(M).value,termSet.INPUT(M).fuzzyset(FormationMatrix(rule,M)).value,3)*Lambda1Set.INPUT(M).fuzzyset(FormationMatrix(rule,M)); %dr/dx
                        temp=r.*exp(j.*(theata1Ofh));
                        testmembership1=testmembership1.*temp;
                        temp2=real(temp);
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
                testNBeta(N).value(1:length(FormationMatrix),1:NumberOfTestPoint)=0;
                for jj=1:NumberOfTestPoint
                    for rule=1:length(FormationMatrix)
                        SumReal=sum(real(testBeta(N).value(:,jj)));
                        SumImag=sum(imag(testBeta(N).value(:,jj)));
                        if SumReal==0
                            SumReal=1e-99;
                        end
                        if SumImag==0
                            SumImag=1e-99;
                        end
                        temp1=real(testBeta(N).value(rule,jj))/SumReal;
                        temp2=imag(testBeta(N).value(rule,jj))/SumImag;
                        testNBeta(N).value(rule,jj)=temp1+temp2*j;
                    end
                end
           end
        
        
        for N=1:NumberOfOUTPUT
            %把每個beta裡的規則分開儲存
            for K=1:NumberOfPremise
                testB(N).k(K).value=testNBeta(N).value(K,:);
            end
        end
        
             
        for jj=1:NumberOfTestPoint      
            temp2=[];
            for N=1:NumberOfOUTPUT
                temp1=[];
                for K=1:NumberOfPremise    
                    testlambda=testNBeta(N).value(K,jj);
                    temp=[];   
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
             testyHead(jj,:)=testA(jj).value*PSOgBest.RLSE.theta;  %y
        end


        %caculate error    
       for N=1:NumberOfOUTPUT
            testError(:,N)=testy(N).value-testyHead(:,N);
       end

            
            