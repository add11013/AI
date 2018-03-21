function time=MODEL_PSO(trial, ExperimentationName, NumberOfTrainPoint, ParameterOfPreSC, ParameterOfConsSC)
close all;
tic

%use function "dataSplit" to split data
OriginalData=DataSplit(ExperimentationName);
NumberOfTestPoint=size(OriginalData,1)-NumberOfTrainPoint;
%use function "FeatureSelection(OriginalData)" to do the multi-target feature
%selection
[FeatureIndex, DataMatrix]=FeatureSelection(OriginalData, NumberOfTrainPoint, ExperimentationName);
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
k=3*NumberOfTarget+1;
for N=1:NumberOfOUTPUT
    realPartOfTrain=DataMatrix(1:NumberOfTrainPoint,k);
    imagPartOfTrain=DataMatrix(1:NumberOfTrainPoint,k+1);
    k=k+2;
    y(N).value=realPartOfTrain+imagPartOfTrain*j;
end

    
%% caculate NumberOfCons
% combine all h, caculate what number of clusters, then use Fuzzy C-mean to
% clustering
TMP0=[];
for M=1:NumberOfINPUT
    TMP=h(M).value;
    TMP0=[TMP0;TMP];
end
TotalINPUT=TMP0;
%Number Of Consequence
NumberOfCons=length(subclust(TotalINPUT,ParameterOfConsSC));

%% Cons. fuzzyset parameters
% get cons. center and standard by using Fuzzy C-mean
for N=1:NumberOfOUTPUT
    ConsCenter(N).value=fcm(y(N).value,NumberOfCons);
    ConsStd(N).value=std(real(y(N).value))+j*std(imag(y(N).value));
    y(N).mean=mean(y(N).value);
    y(N).std=std(real(y(N).value))+std(imag(y(N).value))*j;
end

%% Consequences
%caculate cons. center and standard of each output 
        for N=1:NumberOfOUTPUT
            for Q=1:NumberOfCons
                C(N).q(Q)=ConsCenter(N).value(Q);
                S(N).q(Q)=ConsStd(N).value;
            end
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
upper=15;
down=4;
    if length(PreFormationMatrix)<=upper
        %介於上下界之間，正常
        if length(PreFormationMatrix)>=down
            FormationMatrix=PreFormationMatrix;
            CubeBetaSum=OrderCubeBetaSum;
        %小於下界，拿原始的來補
        else
            FormationMatrix=PreFormationMatrix(1:down,:);
            CubeBetaSum=PreFormationMatrix(1:down);
        end
    %大於上界，只留下前面的
    else
        FormationMatrix=PreFormationMatrix(1:upper,:);
        CubeBetaSum=OrderCubeBetaSum(1:upper);
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
    NumberOfPremiseParameters=NumberOfPremiseParameters+length(h(M).center)*4;
end

%% PSO parameters
  PSO.w=0.8;
  PSO.c1=2;
  PSO.c2=2;
  PSO.s1=rand(1);
  PSO.s2=rand(1);
  PSO.swarm_size=64;
  PSO.iterations=100;
  %initialize the particles
  for i=1:PSO.swarm_size
    j1=1;
    for M=1:NumberOfINPUT
        %每一個FUZZYSET有center、std、lambda1、lambda2
        for ii=1:length(h(M).center)
            particle(i).Position(j1)=h(M).center(ii)+randn; %center
            particle(i).Position(j1+1)=h(M).std+randn; %std
            particle(i).Position(j1+2)=randn; %lambda1
            particle(i).Position(j1+3)=randn; %lambda2
            j1=j1+4;
        end
    end
    particle(i).Velocity(1:NumberOfPremiseParameters)=0;
    particle(i).pBestPosition=particle(i).Position;
    particle(i).pBestDistance=1e99;
  end
  PSOgBest.Position=particle(1).Position;
  PSOgBest.Distance=1e99;
  
%% RLSE parameters
        theata0(1:(NumberOfINPUT+1)*NumberOfCons,1)=0;%
        P0=1e8*eye((NumberOfINPUT+1)*NumberOfCons);
%% main loop
for ite=1:PSO.iterations
  for i=1:PSO.swarm_size
          
      %將粒子的位置，儲存成termSet，以便後續建造前艦部的fuzzy set
            j1=1;
            for M=1:NumberOfINPUT
                for number=1:length(h(M).center)
                    termSet.INPUT(M).fuzzyset(number).value=[particle(i).Position(j1) particle(i).Position(j1+1)];
                    Lambda1Set.INPUT(M).fuzzyset(number)=particle(i).Position(j1+2);
                    Lambda2Set.INPUT(M).fuzzyset(number)=particle(i).Position(j1+3);
                    j1=j1+4;
                end
            end
            
         
         %每個OUTPUT都要算Beta
            %算每一條規則的啟動強度
            for jj=1:NumberOfTrainPoint
                for rule=1:length(FormationMatrix)
                    membership1=1;
                    membership2=1;
%                     membership3=1;
                    for M=1:NumberOfINPUT
                        r=gaussmf(h(M).value(jj),termSet.INPUT(M).fuzzyset(FormationMatrix(rule,M)).value,1);
                        theata1Ofh=gaussmf(h(M).value(jj),termSet.INPUT(M).fuzzyset(FormationMatrix(rule,M)).value,3)*Lambda1Set.INPUT(M).fuzzyset(FormationMatrix(rule,M)); %dr/dx
                        theata2Ofh=gaussmf(h(M).value(jj),termSet.INPUT(M).fuzzyset(FormationMatrix(rule,M)).value,6)*Lambda2Set.INPUT(M).fuzzyset(FormationMatrix(rule,M)); %d2r/dx^2
                        temp=r*exp(j*(theata1Ofh+theata2Ofh));
                        membership1=membership1*temp;
                        temp2=r*cos(theata2Ofh)*cos(theata1Ofh)+r*cos(theata2Ofh)*sin(theata1Ofh)*j;
                        membership2=membership2*temp2;
%                         temp3=r*cos(theata2Ofh)*sin(theata1Ofh)+r*sin(theata2Ofh)*j;
%                         membership3=membership3*temp3;
                    end
                    Beta(1).value(rule,jj)=membership1;
                    Beta(2).value(rule,jj)=membership2;
%                     Beta(3).value(rule,jj)=membership3;
                end
            end
            %Normalization
            for N=1:2
                for jj=1:NumberOfTrainPoint
                    for rule=1:length(FormationMatrix)
                        temp1=real(Beta(N).value(rule,jj))/sum(real(Beta(N).value(:,jj)));
                        temp2=imag(Beta(N).value(rule,jj))/sum(imag(Beta(N).value(:,jj)));
                        NBeta(N).value(rule,jj)=temp1+temp2*j;
                    end
                end
            end

            for N=1:2
                Beta(N).value=NBeta(N).value;
            end
            
        
        
        for N=1:NumberOfOUTPUT
            %把每個beta裡的規則分開儲存
            for K=1:NumberOfPremise
                B(N).k(K).value=Beta(N).value(K,:);
            end
        end
        

        %轉換箭靶中心
        for N=1:NumberOfOUTPUT
            for Q=1:NumberOfCons
                        r2=exp(-(C(N).q(Q)-y(N).mean)*(C(N).q(Q)-y(N).mean)'/(2*y(N).std*(y(N).std)'));
                        w2=r2.*real(-(C(N).q(Q)-y(N).mean)/(y(N).std*(y(N).std)'));
                        aocC(N).q(Q)=r2.*exp(j*w2);
                        r3=exp(-(S(N).q(Q)*(S(N).q(Q))')/(2*y(N).std*(y(N).std)'));
                        w3=r3.*real((-S(N).q(Q))/(y(N).std*(y(N).std)'));
                        aocS(N).q(Q)=r3.*exp(j*w3);
                    for K=1:NumberOfPremise
                        r1=exp(-(B(N).k(K).value-aocC(N).q(Q)).*conj(B(N).k(K).value-aocC(N).q(Q))./(2.*aocS(N).q(Q).*conj(aocS(N).q(Q))));
                        w1=r1.*real(-(B(N).k(K).value-aocC(N).q(Q))/aocS(N).q(Q).*conj(aocS(N).q(Q)));
                        lambda(N).k(K).q(Q).value=r1.*exp(j*w1);
                    end
            end
        end

      %算Q個P
        for Q=1:NumberOfCons
            for jj=1:NumberOfTrainPoint
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
                        TMP=h(M).value(jj);
                        TMP0=[TMP0 TMP];
                    end
                    HOfLambdaTMP=[TMP0];
                    HOfLambdaTMP0=[HOfLambdaTMP0;HOfLambdaTMP];
                end
                HOfLambda=HOfLambdaTMP0;

                P.q(Q).point(jj).value=BOfLambda*LOfLambda*HOfLambda;
            end
        end
        
        

        for k=1:NumberOfTrainPoint
            TMP0=[];
            for Q=1:NumberOfCons
                TMP=P.q(Q).point(k).value;
                TMP0=[TMP0 TMP];
            end
            b.value=transpose(TMP0);
            I=1e-3*eye(NumberOfOUTPUT);
            if(k==1)
                particle(i).RLSE.iteration(k).P=P0-(P0*b.value)/(I+transpose(b.value)*P0*b.value)*transpose(b.value)*P0;
                particle(i).RLSE.iteration(k).theata=theata0+particle(i).RLSE.iteration(k).P*b.value*(y(N).value(k)-transpose(b.value)*theata0);
            else         
                particle(i).RLSE.iteration(k).P=particle(i).RLSE.iteration(k-1).P-(particle(i).RLSE.iteration(k-1).P*b.value)/(I+transpose(b.value)*particle(i).RLSE.iteration(k-1).P*b.value)*transpose(b.value)*particle(i).RLSE.iteration(k-1).P;
                particle(i).RLSE.iteration(k).theata=particle(i).RLSE.iteration(k-1).theata+particle(i).RLSE.iteration(k).P*b.value*(y(N).value(k)-transpose(b.value)*particle(i).RLSE.iteration(k-1).theata);
            end
        end
        for N=1:NumberOfOUTPUT
            TMP0=[];
            TMP=[];
            for Q=1:NumberOfCons
                for jj=1:NumberOfTrainPoint
                    TMP(jj,:)=P.q(Q).point(jj).value(N,:);
                end
                TMP0=[TMP0 TMP];
            end
            particle(i).A(N).value=TMP0;
        end
      %new_yHead(output)
        for N=1:NumberOfOUTPUT
            particle(i).yHead(:,N)=particle(i).A(N).value*particle(i).RLSE.iteration(k).theata;  %y
           %caculate error
            error(:,N)=y(N).value-particle(i).yHead(:,N);
        end

      %累加每一個OUTPUT的RMSE當作RMSE
      temp0=0;
      for N=1:NumberOfOUTPUT
        particle(i).error(:,N)=error(:,N);
        temp=sqrt(sum(error(:,N).*conj(error(:,N)))/(NumberOfTrainPoint));
        temp0=temp0+temp;
      end
      particle(i).rmse=temp0;
      
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
            PSOgBest.RLSE.Theata=particle(i).RLSE.iteration(NumberOfTrainPoint).theata;
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
          plot(1:PSO.iterations,PSO.plotRMSE,'x');
          legend('Learning Curve');
          xlabel('Iteration');
          ylabel('RMSE');

%% test part
 testPart_PSO;
 
%% plot all result              
        %30天漲跌，預測目標是第31天~最後一天的漲跌，也就是最後要判斷第32天~最後一天的實際數值(不是漲跌)
        x=linspace(1,length(OriginalData)-31,length(OriginalData)-31);

        %畫出每個目標
        for i=1:NumberOfTarget
            figure(i)
            hold on
            plot(x,OriginalData(32:length(OriginalData),i));
        end
        
        x=linspace(1,NumberOfAllPoint,NumberOfAllPoint);
        %目標為第31天的漲跌，所以 第31天價格+第31天漲跌(預測出的)=第32天價格(預測出的)
        for N=1:NumberOfOUTPUT
            FinalyHead(:,N)=[PSOgBest.yHead(:,N) ;testyHead(:,N)];
        end
        k=1;

        figure(1)
        output=OriginalData(31:length(OriginalData)-1,1)+real(FinalyHead(:,1));
        plot(x,output);
        line([NumberOfTrainPoint NumberOfTrainPoint],[min(OriginalData(:,1))*0.6 max(OriginalData(:,1))*1.3]);
        legend('Target (TAIEX)','Forecast (TAIEX)');
        xlabel('Trading date index');
        ylabel('Stock price (TAIEX)');
        
        figure(2)
        output=OriginalData(31:length(OriginalData)-1,2)+imag(FinalyHead(:,1));
        plot(x,output);
        line([NumberOfTrainPoint NumberOfTrainPoint],[min(OriginalData(:,2))*0.8 max(OriginalData(:,2))*1.3]);
        legend('Target (DJIA)','Forecast (DJIA)');
        xlabel('Trading date index');
        ylabel('Stock price (DJIA)');
        
        figure(3)
        output=OriginalData(31:length(OriginalData)-1,3)+real(FinalyHead(:,2));
        plot(x,output);
        line([NumberOfTrainPoint NumberOfTrainPoint],[min(OriginalData(:,3))*0.8 max(OriginalData(:,3))*1.3]);
        legend('Target (NASDAQ)','Forecast (NASDAQ)');
        xlabel('Trading date index');
        ylabel('Stock price (NASDAQ)');
        
        figure(4)
        output=OriginalData(31:length(OriginalData)-1,4)+imag(FinalyHead(:,2));
        plot(x,output);
        line([NumberOfTrainPoint NumberOfTrainPoint],[min(OriginalData(:,4))*0.8 max(OriginalData(:,4))*1.3]);
        legend('Target (S&P500)','Forecast (S&P500)');
        xlabel('Trading date index');
        ylabel('Stock price (S&P500)');
        
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
%         RMSETAIEX1=sqrt(sum(real(PSOgBest.Error(:,1)).^2)/(NumberOfTrainPoint));
%         RMSEDJIA1=sqrt(sum(imag(PSOgBest.Error(:,1)).^2)/(NumberOfTrainPoint));
%         RMSETAIEX2=sqrt(sum(real(testError(:,1)).^2)/(NumberOfTestPoint));
%         RMSEDJIA2=sqrt(sum(imag(testError(:,1)).^2)/(NumberOfTestPoint));
%         
%         RMSENASDAQ1=sqrt(sum(real(PSOgBest.Error(:,2)).^2)/(NumberOfTrainPoint));
%         RMSESP1=sqrt(sum(imag(PSOgBest.Error(:,2)).^2)/(NumberOfTrainPoint));
%         RMSENASDAQ2=sqrt(sum(real(testError(:,2)).^2)/(NumberOfTestPoint));
%         RMSESP2=sqrt(sum(imag(testError(:,2)).^2)/(NumberOfTestPoint));
% 
MAPEIBM1=sum(abs(real(PSOgBest.Error(:,1))))/(NumberOfTrainPoint);
MAPEIBM2=sum(abs(real(testError(:,1))))/(NumberOfTestPoint);
MAPEAPPLE1=sum(abs(imag(PSOgBest.Error(:,1))))/(NumberOfTrainPoint);
MAPEAPPLE2=sum(abs(imag(testError(:,1))))/(NumberOfTestPoint);

MAPEDELL1=sum(abs(real(PSOgBest.Error(:,2))))/(NumberOfTrainPoint);
MAPEDELL2=sum(abs(real(testError(:,2))))/(NumberOfTestPoint);
MAPEMicroSoft1=sum(abs(imag(PSOgBest.Error(:,2))))/(NumberOfTrainPoint);
MAPEMicroSoft2=sum(abs(imag(testError(:,2))))/(NumberOfTestPoint);
        
%         RMSE2003_1=sqrt(sum(real(PSOgBest.Error(:,3)).^2)/(NumberOfTrainPoint));
%         RMSE2004_1=sqrt(sum(imag(PSOgBest.Error(:,3)).^2)/(NumberOfTrainPoint));
%         RMSE2003_2=sqrt(sum(real(testError(:,3)).^2)/(NumberOfTestPoint));
%         RMSE2004_2=sqrt(sum(imag(testError(:,3)).^2)/(NumberOfTestPoint));
    time=toc;
    trial=['Result_' ExperimentationName '_trial' int2str(trial)];
    save(trial);
    
 end