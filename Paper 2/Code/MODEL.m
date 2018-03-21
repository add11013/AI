clear
close all;
clc
tic
OriginalData=xlsread('Data_set.csv');
[FeatureIndex DataMatrix]=FeatureSelection(OriginalData);
%% get  Training data
%NumberOfTarget指的是實數型態的目標個數
NumberOfTarget=size(OriginalData,2);
NumberOfTrainPoint=200;
%get h1~hM
for M=1:length(FeatureIndex)
    h(M).value=DataMatrix(1:NumberOfTrainPoint,FeatureIndex(M));
    % substractive clustering for premise fuzzysets
    temp=subclust(h(M).value,0.2);
    h(M).center=sort(temp);
    h(M).std=std(h(M).value);
end


%% prepare target
NumberOfINPUT=length(h);
%NumberOfOUTPUT是指複數型態目標的個數，所以為實數型態除以2
NumberOfOUTPUT=NumberOfTarget/2;

%目標起始行數
k=61;
for N=1:NumberOfOUTPUT
    realPartOfTrain=DataMatrix(1:NumberOfTrainPoint,k);
    imagPartOfTrain=DataMatrix(1:NumberOfTrainPoint,k+1);
    k=k+2;
    y(N).value=realPartOfTrain+imagPartOfTrain*j;
end


%% caculate NumberOfCons
%合併所有h，計算要分幾群之後再用FCM
TMP0=[];
for M=1:NumberOfINPUT
    TMP=h(M).value;
    TMP0=[TMP0;TMP];
end
TotalINPUT=TMP0;
%Number Of Consequence
NumberOfCons=length(subclust(TotalINPUT,0.05));

%% Cons. fuzzyset parameters
%用FCM分群得到後艦部中心和標準差
for N=1:NumberOfOUTPUT
    ConsCenter=fcm(y(N).value,NumberOfCons);
    ConsStd=std(real(y(N).value))+j*std(imag(y(N).value));
    y(N).mean=mean(y(N).value);
    y(N).std=std(y(N).value);
end

%% Consequences
%計算每個OUTPUT的後鑑部中心和標準差
        for N=1:NumberOfOUTPUT
            for Q=1:NumberOfCons
                C(N).q(Q)=ConsCenter(Q);
                S(N).q(Q)=ConsStd;
            end
        end
        
%% formation matrix
%建造一個潛艦部的fuzzyset索引，每一列代表一個cube
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
            formationMatrix(count,M)=k;
            count=count+1;
        end
    end
end

% caculate Number Of Premise Parameters
%累加每個前艦部的參數得到前艦部參數個數
NumberOfPremiseParameters=0;
for M=1:NumberOfINPUT
    %乘以4是因為有center、std、lambda1、lambda2
    NumberOfPremiseParameters=NumberOfPremiseParameters+length(h(M).center)*4;
end

%% firing strength
%將每個INPUT的資料數據帶入算啟動強度
for jj=1:NumberOfTrainPoint
    for cube=1:length(formationMatrix)
        membership=1;
        for M=1:NumberOfINPUT
            temp=gaussmf(h(M).value(jj),[h(M).center(formationMatrix(cube,M)),h(M).std],1);
            membership=membership*temp;    
        end
        BetaOfFormationMatrix(cube,jj)=membership;
    end
end


%% cube selection
for cube=1:length(formationMatrix)
    CubeBetaSum(cube)=sum(BetaOfFormationMatrix(cube,:));
end
CubeBetaMean=mean(CubeBetaSum);
CubeBetaStd=std(CubeBetaSum);
count=1;
%bye為要刪除的cube索引
for cube=1:length(formationMatrix)
    threshold=CubeBetaMean+CubeBetaStd;
    if CubeBetaSum(cube)<threshold
        bye(count)=cube;
        count=count+1;
    end
end
k=0;
for i=1:length(bye)
    formationMatrix(bye(i)-k,:)=[];
    k=k+1;
end
NumberOfPremise=length(formationMatrix);

%% 把formationMatrix沒用到的中心和標準差移除
for M=1:NumberOfINPUT
    k=0;
        for ii=1:length(h(M).center)
            temp=formationMatrix(:,M);
            if ~any(temp==ii)
                h(M).center(ii-k)=[];
                k=k+1;
            end
        end
end

%% 重建一次formationMatrix(重新調整索引)


%% PSO parameters
  PSO.w=0.8;
  PSO.c1=2;
  PSO.c2=2;
  PSO.s1=rand(1);
  PSO.s2=rand(1);
  PSO.swarm_size=100;
  PSO.iterations=10;
  %initialize the particles
  for i=1:PSO.swarm_size
    j1=1;
    for M=1:NumberOfINPUT
        %每一個FUZZYSET有center、std、lambda1、lambda2
        for ii=1:length(h(M).center)
            swarm(i).Position(j1)=h(M).center(ii)+randn*3; %center
            swarm(i).Position(j1+1)=h(M).std+randn*3; %std
            swarm(i).Position(j1+2)=randn; %lambda1
            swarm(i).Position(j1+3)=randn; %lambda2
            j1=j1+4;
        end
    end
    swarm(i).Velocity(1:NumberOfPremiseParameters)=0;
    swarm(i).pBestPosition=swarm(i).Position;
    swarm(i).pBestDistance=1e12;
  end
  PSOgBest.Position=swarm(1).Position;
  PSOgBest.Distance=1e12;
  
  
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
                    termSet.INPUT(M).fuzzyset(number).value=[swarm(i).Position(j1) swarm(i).Position(j1+1)];
                    Lambda1Set.INPUT(M).fuzzyset(number)=swarm(i).Position(j1+2);
                    Lambda2Set.INPUT(M).fuzzyset(number)=swarm(i).Position(j1+3);
                    j1=j1+4;
                end
            end
            
         
         %每個OUTPUT都要算Beta
        for N=1:NumberOfOUTPUT
            Beta(N).value=[];
            %算每一條規則的啟動強度
            for jj=1:NumberOfTrainPoint
                for rule=1:length(formationMatrix)
                    membership=1;
                    for M=1:NumberOfINPUT
                        r=gaussmf(h(M).value(jj),termSet.INPUT(M).fuzzyset(formationMatrix(rule,M)).value,1);
                        theata1Ofh=gaussmf(h(M).value(jj),termSet.INPUT(M).fuzzyset(formationMatrix(rule,M)).value,3)*Lambda1Set.INPUT(M).fuzzyset(formationMatrix(rule,M)); %dr/dx
                        theata2Ofh=gaussmf(h(M).value(jj),termSet.INPUT(M).fuzzyset(formationMatrix(rule,M)).value,6)*Lambda2Set.INPUT(M).fuzzyset(formationMatrix(rule,M)); %d2r/dx^2
                        temp=r*exp(j*(theata1Ofh+theata2Ofh));
                        membership=membership*temp;
                    end
                    Beta(N).value(rule,jj)=membership;
                end
            end
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
                        aocC(N).q(Q)=exp(-(C(N).q(Q)-y(N).mean)*(C(N).q(Q)-y(N).mean)'/(2*y(N).std*(y(N).std)'))  *exp(j*exp(-(C(N).q(Q)-y(N).mean)*(C(N).q(Q)-y(N).mean)'/(2*y(N).std*(y(N).std)'))* real(-(C(N).q(Q)-y(N).mean)/(y(N).std*(y(N).std)')));
                        aocS(N).q(Q)=exp(-(S(N).q(Q)*(S(N).q(Q))')/(2*y(N).std*(y(N).std)'))*exp(j*exp(-(S(N).q(Q)*(S(N).q(Q))')/(2*y(N).std*(y(N).std)'))  *real((-S(N).q(Q))/(y(N).std*(y(N).std)')));
                    for K=1:NumberOfPremise
                         f1(Q).value=exp(-(B(N).k(K).value-aocC(N).q(Q)).*conj(B(N).k(K).value-aocC(N).q(Q))./(2.*aocS(N).q(Q).*conj(aocS(N).q(Q))));
                         lambda(N).k(K).q(Q).value=f1(Q).value.*exp(j.*f1(Q).value.*-(B(N).k(K).value-aocC(N).q(Q))./(((aocS(N).q(Q)).*conj((aocS(N).q(Q))))));
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
                
                P.q(Q).value(jj,:)=[BOfLambda*LOfLambda*HOfLambda];
            end
        end
        
        %
        TMP0=[];
        for Q=1:NumberOfCons
            TMP=P.q(Q).value;
            TMP0=[TMP0 TMP];
        end
        swarm(i).RLSE.A=TMP0;

        b=transpose(swarm(i).RLSE.A);
        for N=1:NumberOfOUTPUT
            for k=1:(NumberOfTrainPoint)
                if(k==1)
                    swarm(i).RLSE.iteration(k).P=P0-(P0*b(:,k)*transpose(b(:,k))*P0)/(1+transpose(b(:,k))*P0*b(:,k));
                    swarm(i).RLSE.iteration(k).theata=theata0+swarm(i).RLSE.iteration(k).P*b(:,k)*(y(N).value(k)-transpose(b(:,k))*theata0);
                else
                    swarm(i).RLSE.iteration(k).P=swarm(i).RLSE.iteration(k-1).P-(swarm(i).RLSE.iteration(k-1).P*b(:,k)*transpose(b(:,k))*swarm(i).RLSE.iteration(k-1).P)/(1+transpose(b(:,k))*swarm(i).RLSE.iteration(k-1).P*b(:,k));
                    swarm(i).RLSE.iteration(k).theata=swarm(i).RLSE.iteration(k-1).theata+swarm(i).RLSE.iteration(k).P*b(:,k)*(y(N).value(k)-transpose(b(:,k))*swarm(i).RLSE.iteration(k-1).theata);
                end
            end
        end
          
      %new_yHead(output)
      for jj=1:NumberOfTrainPoint
        for N=1:NumberOfOUTPUT
            yHead=swarm(i).RLSE.A(jj,:)*swarm(i).RLSE.iteration(k).theata;  %y
            swarm(i).yHead(jj,N)=yHead(N,1);
           %caculate error
            e(jj,N)=(y(N).value(jj)-swarm(i).yHead(jj,N))*conj(y(N).value(jj)-swarm(i).yHead(jj,N));  % target-yHead
            error(jj,N)=y(N).value(jj)-swarm(i).yHead(jj,N);
        end
      end
      %累加每一個OUTPUT的rmse當作rmse
      temp0=0;
      for N=1:NumberOfOUTPUT
        swarm(i).error=error(:,N);
        temp=sqrt(sum(e(:,N))/(NumberOfTrainPoint));
        temp0=temp0+temp;
        swarm(i).rmse=temp0;
      end
      %pbest
        if swarm(i).rmse<swarm(i).pBestDistance
            swarm(i).pBestPosition=swarm(i).Position;        %update pbest position
            swarm(i).pBestDistance=swarm(i).rmse;            %update pbest rmse index
        end
      %gbest
        if swarm(i).rmse<PSOgBest.Distance
            gBest=i;                             %update which one is gbest
            PSOgBest.Distance=swarm(i).rmse;         %update distance of gbest
            PSOgBest.Position=swarm(i).Position;
            PSOgBest.yHead=swarm(i).yHead;
            PSOgBest.Error=swarm(i).error;
            PSOgBest.Theata=swarm(i).RLSE.iteration(NumberOfTrainPoint).theata;
        end

      %update velocity
      swarm(i).Velocity=PSO.w*swarm(i).Velocity+PSO.c1*PSO.s1*(swarm(i).pBestPosition-swarm(i).Position)+PSO.c2*PSO.s2*(PSOgBest.Position-swarm(i).Position);
      %Update swarm Position
      swarm(i).Position=swarm(i).Position+swarm(i).Velocity;  
  end
  PSO.plotRMSE(ite) = PSOgBest.Distance;
  ite
end

%% result
% OUTPUT and Target
      figure(2)
          plot(1:PSO.iterations,PSO.plotRMSE,'x');
          legend('Learning Curve');
          xlabel('Iteration');
          ylabel('RMSE');

%% test part
 testPart
 
%% plot all result      
      
        figure(1)
        hold on
        %30天漲跌，預測目標是第31天~最後一天的漲跌，也就是最後要判斷第32天~最後一天的實際數值(不是漲跌)
        x=linspace(1,length(OriginalData)-31,length(OriginalData)-31);

        %畫出每個目標
        for i=1:NumberOfTarget
            plot(x,OriginalData(32:length(OriginalData),i));
        end
        x=linspace(1,NumberOfAllPoint,NumberOfAllPoint);
        %目標為第31天的漲跌，所以 第31天價格+第31天漲跌(預測出的)=第32天價格(預測出的)
        for N=1:NumberOfOUTPUT
            FinalyHead(:,N)=[PSOgBest.yHead(:,N) ;testyHead(:,N)];
        end
        k=1;
        for N=1:NumberOfOUTPUT
            output=OriginalData(31:291,k)+real(FinalyHead(:,N));
            plot(x,output);
            output=OriginalData(31:291,k+1)+imag(FinalyHead(:,N));
            plot(x,output);
            k=k+2;
        end
        RMSEIBM1=sqrt(sum(real(PSOgBest.Error(:,1)).^2)/(NumberOfTrainPoint));
        RMSETSMC1=sqrt(sum(imag(PSOgBest.Error(:,1)).^2)/(NumberOfTrainPoint));
        RMSEIBM2=sqrt(sum(real(testError(:,1)).^2)/(NumberOfTestPoint));
        RMSETSMC2=sqrt(sum(imag(testError(:,1)).^2)/(NumberOfTestPoint));
        legend('Target (IBM)','Target (TSMC)','Forecast (IBM)','Forecast (TSMC)')
        xlabel('Trading date index');
        ylabel('Stock price (IBM ,TSMC)');

          toc      