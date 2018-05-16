clear
close all;
clc
tic
% OriginalData=xlsread('data_Experimentation3.csv');
OriginalData=E2_dataSplit;
[FeatureIndex DataMatrix]=FeatureSelection(OriginalData);
%% get  Training data
%NumberOfTarget指的是實數型態的目標個數
NumberOfTarget=size(OriginalData,2);
NumberOfTrainPoint=219;
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
k=3*NumberOfTarget+1;
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
NumberOfCons=length(subclust(TotalINPUT,0.1));

%% Cons. fuzzyset parameters
%用FCM分群得到後艦部中心和標準差
for N=1:NumberOfOUTPUT
    ConsCenter(N).value=fcm(y(N).value,NumberOfCons);
    ConsStd(N).value=std(real(y(N).value))+j*std(imag(y(N).value));
    y(N).mean=mean(y(N).value);
    y(N).std=std(real(y(N).value))+std(imag(y(N).value))*j;
end

%% Consequences
%計算每個OUTPUT的後鑑部中心和標準差
        for N=1:NumberOfOUTPUT
            for Q=1:NumberOfCons
                C(N).q(Q)=ConsCenter(N).value(Q);
                S(N).q(Q)=ConsStd(N).value;
            end
        end
        
%% formation matrix
%建造一個前艦部的fuzzyset索引，每一列代表一個cube
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
    for cube=1:length(OriginalFormationMatrix)
        membership=1;
        for M=1:NumberOfINPUT
            temp=gaussmf(h(M).value(jj),[h(M).center(OriginalFormationMatrix(cube,M)),h(M).std],1);
            membership=membership*temp;    
        end
        BetaOfFormationMatrix(cube,jj)=membership;
    end
end

%delete cube with threshold
for cube=1:length(OriginalFormationMatrix)
    OriginalCubeBetaSum(cube)=sum(BetaOfFormationMatrix(cube,:));
end
CubeBetaMean=mean(OriginalCubeBetaSum);
CubeBetaStd=std(OriginalCubeBetaSum);
count=1;
%bye為要刪除的cube索引
for cube=1:length(OriginalFormationMatrix)
    threshold=CubeBetaMean+CubeBetaStd;
    if OriginalCubeBetaSum(cube)<threshold
        bye(count)=cube;
        count=count+1;
    end
end
for i=1:length(OriginalCubeBetaSum)
    CubeIndex(i)=i;
end
    %由大到小排序OriginalCubeBetaSum的順序
    %每回合先把第i個當成最大值，挑出最大值與第i個做交換，做i次即排序完畢
    for i=1:length(CubeIndex)-1
        FPMAX=OriginalCubeBetaSum(i);
        MAXINDEX=i;
        %找出最大值的位置
        for ii=i:length(CubeIndex)
            if OriginalCubeBetaSum(ii)>FPMAX
                FPMAX=OriginalCubeBetaSum(ii);
                MAXINDEX=ii;
            end
        end
        %這回合最大值與第i個p交換
        ptemp=OriginalCubeBetaSum(i);
        OriginalCubeBetaSum(i)=OriginalCubeBetaSum(MAXINDEX);
        OriginalCubeBetaSum(MAXINDEX)=ptemp;

        indextemp=CubeIndex(i);
        CubeIndex(i)=CubeIndex(MAXINDEX);
        CubeIndex(MAXINDEX)=indextemp;
        
        FMtemp=OriginalFormationMatrix(i,:);
        OriginalFormationMatrix(i,:)=OriginalFormationMatrix(MAXINDEX,:);
        OriginalFormationMatrix(MAXINDEX,:)=FMtemp;

    end
PreFormationMatrix=OriginalFormationMatrix;
PreCubeBetaSum=OriginalCubeBetaSum;

k=0;
for i=1:length(bye)
    PreFormationMatrix(bye(i)-k,:)=[];
    PreCubeBetaSum(bye(i)-k)=[];
    k=k+1;
end
upper=15;
down=4;
    if length(PreFormationMatrix)<=upper
        %介於上下界之間，正常
        if length(PreFormationMatrix)>=down
            FormationMatrix=PreFormationMatrix;
            CubeBetaSum=PreCubeBetaSum;
        %小於下界，拿原始的來補
        else
            FormationMatrix=OriginalFormationMatrix(1:down,:);
            CubeBetaSum=OriginalCubeBetaSum(1:down,:);
        end
    %大於上界，只留下前面的
    else
        FormationMatrix=PreFormationMatrix(1:upper,:);
        CubeBetaSum=PreCubeBetaSum(1:upper);
    end
%% 把FormationMatrix沒用到的中心和標準差移除，並且重新編碼
    tag=0;
    for M=1:NumberOfINPUT
        tagM(M)=0;
        k=0;
            for ii=1:length(h(M).center)
                temp=FormationMatrix(:,M);
                %如果不存在FormationMatrix
                if ~any(temp==ii)
                    tag=1;
                    tagM(M)=1;
                    byebye(M,k+1)=ii;
                    %delete center
                    h(M).center(ii-k)=[];
                    k=k+1;
                end
            end
    end

%重新編碼FormationMatrix
    if tag==1
        for M=1:NumberOfINPUT
            if tagM(M)~=0
            for i=1:length(byebye(M,:))
                %
                if byebye(M,i)~=0
                    for ii=1:length(FormationMatrix)
                        if FormationMatrix(ii,M)>byebye(M,i)
                            FormationMatrix(ii,M)=FormationMatrix(ii,M)-1;
                        end
                    end
                end
            end
            end
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
%% RO parameters
  RO.swarm_size=max(3,4*log10(NumberOfPremiseParameters+1));
  RO.iterations=100;
  %initialize the particles
  for i=1:RO.swarm_size
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
  end
  ROgBest.Position=swarm(1).Position;
  ROgBest.Distance=1e24;
  JumpCount=0;
  
%% RLSE parameters
        theata0(1:(NumberOfINPUT+1)*NumberOfCons,1)=0;%
        P0=1e8*eye((NumberOfINPUT+1)*NumberOfCons);

%% main loop
for ite=1:RO.iterations
    stepsize=max(RO.iterations-ite-1.5*NumberOfPremiseParameters,2/NumberOfPremiseParameters)/RO.iterations;
  for i=1:RO.swarm_size
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
            Beta(1).value=[];
            %算每一條規則的啟動強度
            for jj=1:NumberOfTrainPoint
                for rule=1:length(FormationMatrix)
                    membership1=1;
                    membership2=1;
                    for M=1:NumberOfINPUT
                        r=gaussmf(h(M).value(jj),termSet.INPUT(M).fuzzyset(FormationMatrix(rule,M)).value,1);
                        theata1Ofh=gaussmf(h(M).value(jj),termSet.INPUT(M).fuzzyset(FormationMatrix(rule,M)).value,3)*Lambda1Set.INPUT(M).fuzzyset(FormationMatrix(rule,M)); %dr/dx
                        theata2Ofh=gaussmf(h(M).value(jj),termSet.INPUT(M).fuzzyset(FormationMatrix(rule,M)).value,6)*Lambda2Set.INPUT(M).fuzzyset(FormationMatrix(rule,M)); %d2r/dx^2
                        temp=r*exp(j*(theata1Ofh+theata2Ofh));
                        membership1=membership1*temp;
                        temp2=r*sin(theata2Ofh)+r*cos(theata2Ofh)*sin(theata1Ofh)*j;
                        membership2=membership2*temp2;
                    end
                    Beta(1).value(rule,jj)=membership1;
                    Beta(2).value(rule,jj)=membership2;
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
        
        
        for N=1:NumberOfOUTPUT
            TMP0=[];
            TMP=[];
            for Q=1:NumberOfCons
                for jj=1:NumberOfTrainPoint
                    TMP(jj,:)=P.q(Q).point(jj).value(N,:);
                end
                TMP0=[TMP0 TMP];
            end
            swarm(i).RLSE(N).A=TMP0;
        end
        for N=1:NumberOfOUTPUT
            b(N).value=transpose(swarm(i).RLSE(N).A);
        end
        for N=1:NumberOfOUTPUT
            for k=1:(NumberOfTrainPoint)
                if(k==1)
                    swarm(i).RLSE(N).iteration(k).P=P0-(P0*b(N).value(:,k)*transpose(b(N).value(:,k))*P0)/(1+transpose(b(N).value(:,k))*P0*b(N).value(:,k));
                    swarm(i).RLSE(N).iteration(k).theata=theata0+swarm(i).RLSE(N).iteration(k).P*b(N).value(:,k)*(y(N).value(k)-transpose(b(N).value(:,k))*theata0);
                else
                    swarm(i).RLSE(N).iteration(k).P=swarm(i).RLSE(N).iteration(k-1).P-(swarm(i).RLSE(N).iteration(k-1).P*b(N).value(:,k)*transpose(b(N).value(:,k))*swarm(i).RLSE(N).iteration(k-1).P)/(1+transpose(b(N).value(:,k))*swarm(i).RLSE(N).iteration(k-1).P*b(N).value(:,k));
                    swarm(i).RLSE(N).iteration(k).theata=swarm(i).RLSE(N).iteration(k-1).theata+swarm(i).RLSE(N).iteration(k).P*b(N).value(:,k)*(y(N).value(k)-transpose(b(N).value(:,k))*swarm(i).RLSE(N).iteration(k-1).theata);
                end
            end
        end
          
      %new_yHead(output)
        for N=1:NumberOfOUTPUT
            swarm(i).yHead(:,N)=swarm(i).RLSE(N).A*swarm(i).RLSE(N).iteration(k).theata;  %y
           %caculate error
            error(:,N)=y(N).value-swarm(i).yHead(:,N);
        end

      %累加每一個OUTPUT的rmse當作rmse
      temp0=0;
      for N=1:NumberOfOUTPUT
        swarm(i).error(:,N)=error(:,N);
        temp=sqrt(sum(error(:,N).*conj(error(:,N)))/(NumberOfTrainPoint));
        temp0=temp0+temp;
      end
      swarm(i).rmse=temp0;
      %gbest
        if swarm(i).rmse<ROgBest.Distance
            gBest=i;                             %update which one is gbest
            ROgBest.Distance=swarm(i).rmse;         %update distance of gbest
            ROgBest.Position=swarm(i).Position;
            ROgBest.yHead=swarm(i).yHead;
            ROgBest.Error=swarm(i).error;
            for N=1:NumberOfOUTPUT
                ROgBest.RLSE(N).Theata=swarm(i).RLSE(N).iteration(NumberOfTrainPoint).theata;
            end
        else
            JumpCount=JumpCount+1;
        end
        
        if JumpCount>10
            j1=1;
            %每一個FUZZYSET有center、std、lambda1、lambda2
            for ii=1:NumberOfPremiseParameters/4
                swarm(i).Position(j1)=ROgBest.Position(j1)+randn*5*stepsize; %center
                swarm(i).Position(j1+1)=ROgBest.Position(j1+1)+randn*5*stepsize; %std
                swarm(i).Position(j1+2)=ROgBest.Position(j1+2)+randn*stepsize; %lambda1
                swarm(i).Position(j1+3)=ROgBest.Position(j1+3)+randn*stepsize; %lambda2
                j1=j1+4;
            end
            JumpCount=0;
        else
            j1=1;
            %每一個FUZZYSET有center、std、lambda1、lambda2
            for ii=1:NumberOfPremiseParameters/4
                swarm(i).Position(j1)=ROgBest.Position(j1)+randn*3*stepsize; %center
                swarm(i).Position(j1+1)=ROgBest.Position(j1+1)+randn*3*stepsize; %std
                swarm(i).Position(j1+2)=ROgBest.Position(j1+2)+randn*stepsize; %lambda1
                swarm(i).Position(j1+3)=ROgBest.Position(j1+3)+randn*stepsize; %lambda2
                j1=j1+4;
            end
        end


  end
  RO.plotRMSE(ite) = ROgBest.Distance;
  ite
end

%% result
% OUTPUT and Target
      figure(5)
          plot(1:RO.iterations,RO.plotRMSE,'x');
          legend('Learning Curve');
          xlabel('Iteration');
          ylabel('RMSE');

%% test part
 testPart_RO
 
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
            FinalyHead(:,N)=[ROgBest.yHead(:,N) ;testyHead(:,N)];
        end
        k=1;

        figure(1)
        output=OriginalData(31:282,1)+real(FinalyHead(:,1));
        plot(x,output);
        legend('Target (IBM)','Forecast (IBM)');
        xlabel('Trading date index');
        ylabel('Stock price (IBM)');
        
        figure(2)
        output=OriginalData(31:282,2)+imag(FinalyHead(:,1));
        plot(x,output);
        legend('Target (TSMC)','Forecast (TSMC)');
        xlabel('Trading date index');
        ylabel('Stock price (TSMC)');
        
        figure(3)
        output=OriginalData(31:282,3)+real(FinalyHead(:,2));
        plot(x,output);
        legend('Target (S&P500)','Forecast (S&P500)');
        xlabel('Trading date index');
        ylabel('Stock price (S&P500)');
        
        figure(4)
        output=OriginalData(31:282,4)+imag(FinalyHead(:,2));
        plot(x,output);
        legend('Target (DJIA)','Forecast (DJIA)');
        xlabel('Trading date index');
        ylabel('Stock price (DJIA)');
        
        RMSE2001_1=sqrt(sum(real(ROgBest.Error(:,1)).^2)/(NumberOfTrainPoint));
        RMSE2002_1=sqrt(sum(imag(ROgBest.Error(:,1)).^2)/(NumberOfTrainPoint));
        RMSE2001_2=sqrt(sum(real(testError(:,1)).^2)/(NumberOfTestPoint));
        RMSE2002_2=sqrt(sum(imag(testError(:,1)).^2)/(NumberOfTestPoint));
        
        RMSE2003_1=sqrt(sum(real(ROgBest.Error(:,2)).^2)/(NumberOfTrainPoint));
        RMSE2004_1=sqrt(sum(imag(ROgBest.Error(:,2)).^2)/(NumberOfTrainPoint));
        RMSE2003_2=sqrt(sum(real(testError(:,2)).^2)/(NumberOfTestPoint));
        RMSE2004_2=sqrt(sum(imag(testError(:,2)).^2)/(NumberOfTestPoint));



          toc      