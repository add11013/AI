function time=MODEL_PSO(trial, ExperimentationName, NumberOfTrainPoint, ParameterOfPreSC)
close all;
tic

OriginalData=DataSplit(ExperimentationName);
%use function "FeatureSelection(OriginalData)" to do the multi-target feature
%selection
[FeatureIndex, DataMatrix]=FeatureSelection(OriginalData, NumberOfTrainPoint, ExperimentationName);

NumberOfAllPoint=size(DataMatrix,1);
NumberOfTestPoint=NumberOfAllPoint-NumberOfTrainPoint;

%% get  Training data
%NumberOfTarget is number of real-value targets
NumberOfTarget=size(OriginalData,2);
%get h1~hM
for M=1:length(FeatureIndex)
    h(M).value=DataMatrix(1:NumberOfTrainPoint,FeatureIndex(M));
    % substractive clustering for premise fuzzysets
    temp=subclust(h(M).value,ParameterOfPreSC);
    h(M).center=sort(temp);
    h(M).std=std(h(M).value);
end


%% prepare target
NumberOfINPUT=length(h);
%NumberOfOUTPUT is number of complex-value targets
NumberOfOUTPUT=NumberOfTarget/2;

%k is the first column of target in DataMatrix 
k=30*NumberOfTarget+1;
if NumberOfTarget==1
    y(1).value=DataMatrix(1:NumberOfTrainPoint,k);
elseif mod(NumberOfTarget,2)==0
    for N=1:NumberOfOUTPUT
        realPartOfTrain=DataMatrix(1:NumberOfTrainPoint,k);
        imagPartOfTrain=DataMatrix(1:NumberOfTrainPoint,k+1);
        k=k+2;
        y(N).value=realPartOfTrain+imagPartOfTrain*j;
    end
else
    for N=1:NumberOfOUTPUT
        realPartOfTrain=DataMatrix(1:NumberOfTrainPoint,k);
        imagPartOfTrain=DataMatrix(1:NumberOfTrainPoint,k+1);
        k=k+2;
        y(N).value=realPartOfTrain+imagPartOfTrain*j;
    end
    y(N+1).value=DataMatrix(1:NumberOfTrainPoint,k);
end
        
%% formation matrix
%construct a premise fuzzyset index, every row is a cube
NumberOfCube=1;
for M=1:NumberOfINPUT
    NumberOfCube=NumberOfCube*length(h(M).center);
end
xxx=1;
yyy=NumberOfCube;
for M=1:NumberOfINPUT
    count=1;
    xxx=xxx*length(h(M).center);
    yyy=yyy/length(h(M).center);
    for i=1:xxx
        for ii=1:yyy
            if mod(i,length(h(M).center))==0
                k=length(h(M).center);
            else
                k=mod(i,length(h(M).center));
            end
            OriginalFormationMatrix(count,M)=k;
            count=count+1;
        end
    end
end


%% cube selection
% firing strength
%將每個INPUT的資料數據帶入算啟動強度
for jj=1:NumberOfTrainPoint
    for cube=1:size(OriginalFormationMatrix,1)
        membership=1;
        for M=1:NumberOfINPUT
            temp=gaussmf(h(M).value(jj),[h(M).center(OriginalFormationMatrix(cube,M)),h(M).std],1);
            membership=membership*temp;
        end
        BetaOfFormationMatrix(cube,jj)=membership;
    end
end
for cube=1:size(OriginalFormationMatrix,1)
    OriginalCubeBetaSum(cube)=sum(BetaOfFormationMatrix(cube,:));
end

OrderFormationMatrix=OriginalFormationMatrix;
OrderCubeBetaSum=OriginalCubeBetaSum;
    %由大到小排序OriginalCubeBetaSum的順序
    for i=1:length(OrderFormationMatrix)
        for ii=1:length(OrderFormationMatrix)
            if OrderCubeBetaSum(ii)<OrderCubeBetaSum(i)
                    ptemp=OrderCubeBetaSum(i);
                    OrderCubeBetaSum(i)=OrderCubeBetaSum(ii);
                    OrderCubeBetaSum(ii)=ptemp;

                    FMtemp=OrderFormationMatrix(i,:);
                    OrderFormationMatrix(i,:)=OrderFormationMatrix(ii,:);
                    OrderFormationMatrix(ii,:)=FMtemp;
            end
        end
    end
    
%decide which cube should be deleted by threshold
CubeBetaMean=mean(OrderCubeBetaSum);
CubeBetaStd=std(OrderCubeBetaSum);
count=1;
%bye為要刪除的cube索引
for cube=1:size(OriginalFormationMatrix,1)
    threshold=CubeBetaMean+CubeBetaStd;
    if OrderCubeBetaSum(cube)<threshold
        bye(count)=cube;
        count=count+1;
    end
end


k=1;
for i=1:length(OriginalFormationMatrix)
    same=0;
    for ii=1:length(bye)
        if i==bye(ii)
            same=1;
        end
    end
    if same==0
        PreFormationMatrix(k,:)=OrderFormationMatrix(i,:);
        k=k+1;
    end
end
RuleUpper=15;
RuleDown=4;
    if length(PreFormationMatrix)<=RuleUpper
        %介於上下界之間，正常
        if length(PreFormationMatrix)>=RuleDown
            FormationMatrix=PreFormationMatrix;
            CubeBetaSum=OrderCubeBetaSum;
        %小於下界，拿原始的來補
        else
            FormationMatrix=PreFormationMatrix(1:RuleDown,:);
            CubeBetaSum=PreFormationMatrix(1:RuleDown);
        end
    %大於上界，只留下前面的
    else
        FormationMatrix=PreFormationMatrix(1:RuleUpper,:);
        CubeBetaSum=OrderCubeBetaSum(1:RuleUpper);
    end
%% 把FormationMatrix沒用到的中心和標準差移除，並且重新編碼
    for M=1:NumberOfINPUT
        k=0;
            for ii=1:length(h(M).center)
                temp=FormationMatrix(:,M);
                if ~any(temp==ii)
                    %delete center
                    h(M).center(ii-k)=[];
                    k=k+1;
                end
            end
    end

%重新編碼formationMatrix，建造一個前艦部的fuzzyset索引，每一列代表一個cube
for M=1:NumberOfINPUT
    FMmin=min(FormationMatrix(:,M));
    for ii=1:length(FormationMatrix)
        FormationMatrix(ii,M)=FormationMatrix(ii,M)-FMmin+1;
    end
end

    
NumberOfPremise=length(FormationMatrix);
% caculate Number Of Premise Parameters
%累加每個前艦部的參數得到前艦部參數個數
NumberOfPremiseParameters=0;
for M=1:NumberOfINPUT
    %乘以4是因為有center、std、lambda1、lambda2
    NumberOfPremiseParameters=NumberOfPremiseParameters+length(h(M).center)*3;
end
toc
%% PSO parameters
  PSO.w=0.8;
  PSO.c1=2;
  PSO.c2=2;
  PSO.s1=rand(1);
  PSO.s2=rand(1);
  PSO.swarm_size=50;
  PSO.iterations=100;
  %initialize the particles
  for i=1:PSO.swarm_size
    j1=1;
    for M=1:NumberOfINPUT
        %每一個FUZZYSET有center、std、lambda1、lambda2
        for ii=1:length(h(M).center)
            particle(i).Position(j1)=h(M).center(ii)^2*randn; %center
            particle(i).Position(j1+1)=h(M).std^2*randn; %std
            particle(i).Position(j1+2)=randn*10; %lambda1
            j1=j1+3;
        end
    end
    particle(i).Velocity(1:NumberOfPremiseParameters)=0;
    particle(i).pBestPosition=particle(i).Position;
    particle(i).pBestDistance=1e99;
    particle(i).yHead(1:NumberOfTrainPoint,1:NumberOfOUTPUT)=0;
    particle(i).error(1:NumberOfTrainPoint,1:NumberOfOUTPUT)=1e99;
  end
  PSOgBest.Position=particle(1).Position;
  PSOgBest.Distance=1e99;
  
%% RLSE parameters
        NumberOfCons=NumberOfPremise;
        theta0(1:(NumberOfINPUT+1)*NumberOfCons,1)=0;%
        P0=1e8*eye((NumberOfINPUT+1)*NumberOfCons);
%% main loop
tic
for ite=1:PSO.iterations
  for i=1:PSO.swarm_size
          
      %將粒子的位置，儲存成termSet，以便後續建造前艦部的fuzzy set
            j1=1;
            for M=1:NumberOfINPUT
                for number=1:length(h(M).center)
                    termSet.INPUT(M).fuzzyset(number).value=[particle(i).Position(j1) particle(i).Position(j1+1)];
                    Lambda1Set.INPUT(M).fuzzyset(number)=particle(i).Position(j1+2);
                    j1=j1+3;
                end
            end

         
         %每個OUTPUT都要算Beta
            %算每一條規則的啟動強度

                for rule=1:length(FormationMatrix)
                    membership1=1;
                    membership2=1;
%                     membership3=1;
                    for M=1:NumberOfINPUT
                        r=gaussmf(h(M).value,termSet.INPUT(M).fuzzyset(FormationMatrix(rule,M)).value,1);
                        theta1Ofh=gaussmf(h(M).value,termSet.INPUT(M).fuzzyset(FormationMatrix(rule,M)).value,3)*Lambda1Set.INPUT(M).fuzzyset(FormationMatrix(rule,M)); %dr/dx
                        temp=r.*exp(j.*(theta1Ofh));
                        membership1=membership1.*temp;
                        temp2=real(temp);
                        membership2=membership2.*temp2;
%                         temp3=r*cos(theta2Ofh)*sin(theta1Ofh)+r*sin(theta2Ofh)*j;
%                         membership3=membership3*temp3;
                    end
                    Beta(1).value(rule,:)=membership1;
                    Beta(2).value(rule,:)=membership2;
%                     Beta(3).value(rule,jj)=membership3;
                end


            %Normalization
            for N=1:NumberOfOUTPUT
                NBeta(N).value(1:length(FormationMatrix),1:NumberOfTrainPoint)=0;
                for jj=1:NumberOfTrainPoint
                    for rule=1:length(FormationMatrix)
                        SumReal=sum(real(Beta(N).value(:,jj)));
                        SumImag=sum(imag(Beta(N).value(:,jj)));
                        if SumReal==0
                            SumReal=1e-99;
                        end
                        if SumImag==0
                            SumImag=1e-99;
                        end
                        temp1=real(Beta(N).value(rule,jj))/SumReal;
                        temp2=imag(Beta(N).value(rule,jj))/SumImag;
                        NBeta(N).value(rule,jj)=temp1+temp2*j;
                    end
                end
            end

        
        for N=1:NumberOfOUTPUT
            %把每個beta裡的規則分開儲存
            for K=1:NumberOfPremise
                B(N).k(K).value=NBeta(N).value(K,:);
            end
        end
        
           
        for jj=1:NumberOfTrainPoint      
            temp2=[];
            for N=1:NumberOfOUTPUT
                temp1=[];
                for K=1:NumberOfPremise    
                    lambda=NBeta(N).value(K,jj);
                    temp=[];
                    temp=lambda;
                    for M=1:NumberOfINPUT
                        temp=[temp lambda*h(M).value(jj)];
                    end
                    temp1=[temp1 temp];
                end
                temp2=[temp2; temp1];                
            end
            b(jj).value=transpose(temp2);
        end

 
        %RLSE
        for k=1:NumberOfTrainPoint
            I=1e-9*eye(NumberOfOUTPUT);
            if(k==1)
                particle(i).RLSE.P=P0-(P0*b(k).value)/(I+transpose(b(k).value)*P0*b(k).value)*transpose(b(k).value)*P0;
                particle(i).RLSE.theta=theta0+particle(i).RLSE.P*b(k).value*(y(N).value(k)-transpose(b(k).value)*theta0);
            else
                particle(i).RLSE.P=particle(i).RLSE.P-(particle(i).RLSE.P*b(k).value)/(I+transpose(b(k).value)*particle(i).RLSE.P*b(k).value)*transpose(b(k).value)*particle(i).RLSE.P;
                particle(i).RLSE.theta=particle(i).RLSE.theta+particle(i).RLSE.P*b(k).value*(y(N).value(k)-transpose(b(k).value)*particle(i).RLSE.theta);
            end
        end
            
      %new_yHead(output)
        for jj=1:NumberOfTrainPoint
             particle(i).yHead(jj,:)=transpose(b(jj).value)*particle(i).RLSE.theta;  %y
        end
       %caculate error
       for N=1:NumberOfOUTPUT
            error(:,N)=y(N).value-particle(i).yHead(:,N);
       end
      %累加每一個OUTPUT的RMSE當作RMSE
      temp0=0;
      particle(i).error=error;
      for N=1:NumberOfOUTPUT
        temp=sqrt(sum(error(:,N).*conj(error(:,N)))/(NumberOfTrainPoint));
        temp0=temp0+temp;
      end
      particle(i).rmse=temp0/NumberOfOUTPUT;
      
        %pbest
        if particle(i).rmse<particle(i).pBestDistance
            particle(i).pBestPosition=particle(i).Position;        %update pbest position
            particle(i).pBestDistance=particle(i).rmse;            %update pbest rmse index
        end
      %gbest
        if particle(i).rmse<PSOgBest.Distance
            gBest=i;                             %update which one is gbest
            PSOgBest.Distance=particle(i).rmse;         %update distance of gbest
            PSOgBest.Position=particle(i).Position;
            PSOgBest.yHead=particle(i).yHead;
            PSOgBest.Error=particle(i).error;
            PSOgBest.RLSE.theta=particle(i).RLSE.theta;
        end
        
      %update velocity
      particle(i).Velocity=PSO.w*particle(i).Velocity+PSO.c1*PSO.s1*(particle(i).pBestPosition-particle(i).Position)+PSO.c2*PSO.s2*(PSOgBest.Position-particle(i).Position);
      %Update particle Position
      particle(i).Position=particle(i).Position+particle(i).Velocity;  

  end
  PSO.plotRMSE(ite) = PSOgBest.Distance;
  ite

end

%% result
% OUTPUT and Target
          figure(7)
          plot(1:PSO.iterations,PSO.plotRMSE);
          legend('Learning Curve');
          xlabel('Iteration');
          ylabel('RMSE');

%% test part
 testPart_PSO;
 
% plot all result              
%        30天漲跌，預測目標是第31天~最後一天的漲跌，也就是最後要判斷第32天~最後一天的實際數值(不是漲跌)
        x=linspace(1,length(OriginalData)-31,length(OriginalData)-31);

%         畫出每個目標
        for i=1:NumberOfTarget
            figure(i)
            hold on
            plot(x,OriginalData(32:length(OriginalData),i));
        end
        
        x=linspace(2,NumberOfAllPoint+1,NumberOfAllPoint);
%         目標為第31天的漲跌，所以 第31天價格+第31天漲跌(預測出的)=第32天價格(預測出的)
        for N=1:NumberOfOUTPUT
            FinalyHead(:,N)=[PSOgBest.yHead(:,N) ;testyHead(:,N)];
        end
        k=1;

        figure(1)
        output=OriginalData(31:length(OriginalData)-1,1)+real(FinalyHead(:,1));
        plot(x,output);
        line([NumberOfTrainPoint NumberOfTrainPoint],[min(OriginalData(:,1))*0.6 max(OriginalData(:,1))*1.3]);
        legend('Target (APPLE)','Forecast (APPLE)');
        xlabel('Trading date index');
        ylabel('Stock price (APPLE)');
        
        figure(2)
        output=OriginalData(31:length(OriginalData)-1,2)+imag(FinalyHead(:,1));
        plot(x,output);
        line([NumberOfTrainPoint NumberOfTrainPoint],[min(OriginalData(:,2))*0.8 max(OriginalData(:,2))*1.3]);
        legend('Target (IBM)','Forecast (IBM)');
        xlabel('Trading date index');
        ylabel('Stock price (IBM)');
        
%         figure(3)
%         output=OriginalData(31:length(OriginalData)-1,3)+real(FinalyHead(:,2));
%         plot(x,output);
%         line([NumberOfTrainPoint NumberOfTrainPoint],[min(OriginalData(:,3))*0.8 max(OriginalData(:,3))*1.3]);
%         legend('Target (S&P500)','Forecast (S&P500)');
%         xlabel('Trading date index');
%         ylabel('Stock price (S&P500)');
%         
%         figure(4)
%         output=OriginalData(31:length(OriginalData)-1,4)+imag(FinalyHead(:,2));
%         plot(x,output);
%         line([NumberOfTrainPoint NumberOfTrainPoint],[min(OriginalData(:,4))*0.8 max(OriginalData(:,4))*1.3]);
%         legend('Target (RUSSELL 2000)','Forecast (RUSSELL 2000)');
%         xlabel('Trading date index');
%         ylabel('Stock price (RUSSELL 2000)');
        
%         figure(5)
%         output=OriginalData(31:283,5)+imag(FinalyHead(:,3));
%         plot(x,output);
%         legend('Target (DJIA)','Forecast (DJIA)');
%         xlabel('Trading date index');
%         ylabel('Stock price (DJIA)');
%         
%         figure(6)
%         output=OriginalData(31:283,6)+imag(FinalyHead(:,3));
%         plot(x,output);
%         legend('Target (DJIA)','Forecast (DJIA)');
%         xlabel('Trading date index');
%         ylabel('Stock price (DJIA)');
%         
        
        TrainRMSE1=sqrt(sum(real(PSOgBest.Error(:,1)).^2)/(NumberOfTrainPoint));
        TrainRMSE2=sqrt(sum(imag(PSOgBest.Error(:,1)).^2)/(NumberOfTrainPoint));
        TestRMSE1=sqrt(sum(real(testError(:,1)).^2)/(NumberOfTestPoint));
        TestRMSE2=sqrt(sum(imag(testError(:,1)).^2)/(NumberOfTestPoint));
       
%         TrainRMSE3=sqrt(sum(real(PSOgBest.Error(:,2)).^2)/(NumberOfTrainPoint));
%         TrainRMSE4=sqrt(sum(imag(PSOgBest.Error(:,2)).^2)/(NumberOfTrainPoint));
%         TestRMSE3=sqrt(sum(real(testError(:,2)).^2)/(NumberOfTestPoint));
%         TestRMSE4=sqrt(sum(imag(testError(:,2)).^2)/(NumberOfTestPoint));

% TrainMAPE1=100/NumberOfTrainPoint*sum(abs(real(PSOgBest.Error(:,1)))./real(y(1).value(1:NumberOfTrainPoint)));
% TestMAPE1=100/NumberOfTestPoint*sum(abs(real(testError(:,1)))./real(y(1).value(1:NumberOfTestPoint)));
% TrainMAPE2=100/NumberOfTrainPoint*sum(abs(imag(PSOgBest.Error(:,1)))./(imag(y(1).value(1:NumberOfTrainPoint))));
% TestMAPE2=100/NumberOfTestPoint*sum(abs(imag(testError(:,1)))./imag(y(1).value(1:NumberOfTestPoint)));

% TrainMAPE3=100/NumberOfTrainPoint*sum(abs(real(PSOgBest.Error(:,2)))./real(y(2).value(1:NumberOfTrainPoint)));
% TestMAPE3=100/NumberOfTestPoint*sum(abs(real(testError(:,2)))./real(y(2).value(1:NumberOfTestPoint)));
% TrainMAPE4=100/NumberOfTrainPoint*sum(abs(imag(PSOgBest.Error(:,2)))./(imag(y(2).value(1:NumberOfTrainPoint))));
% TestMAPE4=100/NumberOfTestPoint*sum(abs(imag(testError(:,2)))./imag(y(2).value(1:NumberOfTestPoint)));
         
    time=toc;
    trial=['PSOResult_' ExperimentationName '_trial' int2str(trial)];
    save(trial);
    
 end