function time=MODEL_ACO(trial, ExperimentationName, NumberOfTrainPoint, ParameterOfPreSC, ParameterOfConsSC)
close all;
tic
load(['Data_' ExperimentationName]);

%use function "FeatureSelection(OriginalData)" to do the multi-target feature
%selection
[FeatureIndex, DataMatrix]=FeatureSelection(OriginalData, NumberOfTrainPoint, ExperimentationName);
NumberOfTestPoint=size(DataMatrix,1)-NumberOfTrainPoint;
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
for N=1:NumberOfOUTPUT
    realPartOfTrain=DataMatrix(1:NumberOfTrainPoint,k);
    imagPartOfTrain=DataMatrix(1:NumberOfTrainPoint,k+1);
    k=k+2;
    y(N).value=realPartOfTrain+imagPartOfTrain*j;
end

    
%% calculate NumberOfCons
% combine all h, calculate what number of clusters, then use Fuzzy C-mean to
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
%calculate cons. center and standard of each output 
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
    threshold=CubeBetaMean+1.5*CubeBetaStd;
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
% calculate Number Of Premise Parameters
%累加每個前艦部的參數得到前艦部參數個數
NumberOfPremiseParameters=0;
for M=1:NumberOfINPUT
    %乘以4是因為有center、std、lambda1、lambda2
    NumberOfPremiseParameters=NumberOfPremiseParameters+length(h(M).center)*4;
end
toc
%% ACO parameters
  ACO.q=0.01;
  ACO.convergence=0.85;
  ACO.swarm_size=64;
  ACO.iterations=100;


  ACOgBest.rmse=1e9;
  %initialize the ants
  for i=1:ACO.swarm_size
    j1=1;
    for M=1:NumberOfINPUT
        %每一個FUZZYSET有center、std、lambda1、lambda2
        for ii=1:length(h(M).center)
            ant(i).Position(j1)=h(M).center(ii)+randn; %center
            ant(i).Position(j1+1)=h(M).std+randn; %std
            ant(i).Position(j1+2)=randn; %lambda1
            ant(i).Position(j1+3)=randn; %lambda2
            j1=j1+4;
        end
    end
    ant(i).yHead(1:NumberOfTrainPoint,1:NumberOfOUTPUT)=0;
    ant(i).error(1:NumberOfTrainPoint,1:NumberOfOUTPUT)=1e99;
  end
 
  
%% RLSE parameters
        NumberOfCons=NumberOfPremise;
        theta0(1:(NumberOfINPUT+1)*NumberOfCons,1)=0;%
        P0=1e8*eye((NumberOfINPUT+1)*NumberOfCons);
%% main loop
tic
for ite=1:ACO.iterations
  for i=1:ACO.swarm_size
          
      %將粒子的位置，儲存成termSet，以便後續建造前艦部的fuzzy set
            j1=1;
            for M=1:NumberOfINPUT
                for number=1:length(h(M).center)
                    termSet.INPUT(M).fuzzyset(number).value=[ant(i).Position(j1) ant(i).Position(j1+1)];
                    Lambda1Set.INPUT(M).fuzzyset(number)=ant(i).Position(j1+2);
                    Lambda2Set.INPUT(M).fuzzyset(number)=ant(i).Position(j1+3);
                    j1=j1+4;
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
                        theta2Ofh=gaussmf(h(M).value,termSet.INPUT(M).fuzzyset(FormationMatrix(rule,M)).value,6)*Lambda2Set.INPUT(M).fuzzyset(FormationMatrix(rule,M)); %d2r/dx^2
                        temp=r.*exp(j.*(theta1Ofh+theta2Ofh));
                        membership1=membership1.*temp;
                        temp2=r.*cos(theta2Ofh).*cos(theta1Ofh)+r.*cos(theta2Ofh).*sin(theta1Ofh).*j;
                        membership2=membership2.*temp2;
%                         temp3=r*cos(theta2Ofh)*sin(theta1Ofh)+r*sin(theta2Ofh)*j;
%                         membership3=membership3*temp3;
                    end
                    Beta(1).value(rule,:)=membership1;
                    Beta(2).value(rule,:)=membership2;
%                     Beta(3).value(rule,jj)=membership3;
                end


            %Normalization
            for N=1:2
                NBeta(N).value(1:length(FormationMatrix),1:NumberOfTrainPoint)=0;
                for jj=1:NumberOfTrainPoint
                    for rule=1:length(FormationMatrix)
                        temp1=real(Beta(N).value(rule,jj))/sum(real(Beta(N).value(:,jj)));
                        temp2=imag(Beta(N).value(rule,jj))/sum(imag(Beta(N).value(:,jj)));
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
        
        temp=[];        
        for jj=1:NumberOfTrainPoint      
            temp2=[];
            for N=1:NumberOfOUTPUT
                temp1=[];
                for K=1:NumberOfPremise    
                    lambda=NBeta(N).value(K,jj);
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
            I=1e-5*eye(NumberOfOUTPUT);
            if(k==1)
                ant(i).RLSE.P=P0-(P0*b(k).value)/(I+transpose(b(k).value)*P0*b(k).value)*transpose(b(k).value)*P0;
                ant(i).RLSE.theta=theta0+ant(i).RLSE.P*b(k).value*(y(N).value(k)-transpose(b(k).value)*theta0);
            else         
                ant(i).RLSE.P=ant(i).RLSE.P-(ant(i).RLSE.P*b(k).value)/(I+transpose(b(k).value)*ant(i).RLSE.P*b(k).value)*transpose(b(k).value)*ant(i).RLSE.P;
                ant(i).RLSE.theta=ant(i).RLSE.theta+ant(i).RLSE.P*b(k).value*(y(N).value(k)-transpose(b(k).value)*ant(i).RLSE.theta);
            end
        end

            A=b;
            
      %new_yHead(output)
        for jj=1:NumberOfTrainPoint
             ant(i).yHead(jj,:)=transpose(A(jj).value)*ant(i).RLSE.theta;  %y
        end
       %calculate error
       for N=1:NumberOfOUTPUT
            error(:,N)=y(N).value-ant(i).yHead(:,N);
       end
      %累加每一個OUTPUT的RMSE當作RMSE
      temp0=0;
      ant(i).error=error;
      for N=1:NumberOfOUTPUT
        temp=sqrt(sum(error(:,N).*conj(error(:,N)))/(NumberOfTrainPoint));
        temp0=temp0+temp;
      end
      ant(i).rmse=temp0/NumberOfOUTPUT;
  end
    SortMatrix(1:ACO.swarm_size,1:5)=0;
  %SortMatrix(:,1)=>RMSE
  %SortMatrix(:,2)=>no. of ant
  %SortMatrix(:,3)=>sequence
  %SortMatrix(:,4)=>weight
  %SortMatrix(:,5)=>probability
  TMP=linspace(1,ACO.swarm_size,ACO.swarm_size);
  SortMatrix=arrayfun(@(x) x.rmse,ant);
  SortMatrix=transpose([SortMatrix;TMP]);
  SortMatrix=sortrows(SortMatrix,1);
  SortMatrix=[SortMatrix transpose(TMP)];
  ACO.l=SortMatrix(:,3);
  SortMatrix(:,4)=CalculateWeight(ACO.q, ACO.swarm_size,ACO.l);
  SortMatrix(:,5)=SortMatrix(:,4)./sum(SortMatrix(:,4));
  Selected=Roulette(SortMatrix(:,5));
  SelectedAnt=SortMatrix(Selected,2);
  
  TMP=[];
  TMP(1:ACO.swarm_size,1:NumberOfPremiseParameters)=0;
  for i=1:ACO.swarm_size
      TMP(i,:)=ant(i).Position;
  end
  AntMean=sum(TMP)./ACO.swarm_size;
  
  for cc=1:NumberOfPremiseParameters
          temp=0;
          for ii=1:ACO.swarm_size
              temp=temp+abs(TMP(ii,cc)-TMP(SelectedAnt,cc))/(ACO.swarm_size-1);
          end
          SelectedStd(cc)=temp;
  end
  
  for i=(ACO.swarm_size/2+1):ACO.swarm_size
      replaced=SortMatrix(i,2);
            ant(replaced).Position=randn*SelectedStd+ant(SelectedAnt).Position;
  end
  
  if ant(SelectedAnt).rmse<ACOgBest.rmse
      ACOgBest.rmse=ant(SelectedAnt).rmse;
      ACOgBest.Position=ant(SelectedAnt).Position;
      ACOgBest.yHead=ant(SelectedAnt).yHead;
      ACOgBest.Error=ant(SelectedAnt).error;
      ACOgBest.theta=ant(i).RLSE.theta;
      
  end
  ACO.plotRMSE(ite) = ACOgBest.rmse;
  ite

end

%% result
% OUTPUT and Target
      figure(7)
          plot(1:ACO.iterations,ACO.plotRMSE);
          legend('Learning Curve');
          xlabel('Iteration');
          ylabel('RMSE');

%% test part
 testPart_ACO;
 
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
            FinalyHead(:,N)=[ACOgBest.yHead(:,N) ;testyHead(:,N)];
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
        
        figure(3)
        output=OriginalData(31:length(OriginalData)-1,3)+real(FinalyHead(:,2));
        plot(x,output);
        line([NumberOfTrainPoint NumberOfTrainPoint],[min(OriginalData(:,3))*0.8 max(OriginalData(:,3))*1.3]);
        legend('Target (S&P500)','Forecast (S&P500)');
        xlabel('Trading date index');
        ylabel('Stock price (S&P500)');
        
        figure(4)
        output=OriginalData(31:length(OriginalData)-1,4)+imag(FinalyHead(:,2));
        plot(x,output);
        line([NumberOfTrainPoint NumberOfTrainPoint],[min(OriginalData(:,4))*0.8 max(OriginalData(:,4))*1.3]);
        legend('Target (RUSSELL 2000)','Forecast (RUSSELL 2000)');
        xlabel('Trading date index');
        ylabel('Stock price (RUSSELL 2000)');
        
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
        TrainRMSE1=sqrt(sum(real(ACOgBest.Error(:,1)).^2)/(NumberOfTrainPoint));
        TrainRMSE2=sqrt(sum(imag(ACOgBest.Error(:,1)).^2)/(NumberOfTrainPoint));
        TestRMSE1=sqrt(sum(real(testError(:,1)).^2)/(NumberOfTestPoint));
        TestRMSE2=sqrt(sum(imag(testError(:,1)).^2)/(NumberOfTestPoint));
        
        TrainRMSE3=sqrt(sum(real(ACOgBest.Error(:,2)).^2)/(NumberOfTrainPoint));
        TrainRMSE4=sqrt(sum(imag(ACOgBest.Error(:,2)).^2)/(NumberOfTrainPoint));
        TestRMSE3=sqrt(sum(real(testError(:,2)).^2)/(NumberOfTestPoint));
        TestRMSE4=sqrt(sum(imag(testError(:,2)).^2)/(NumberOfTestPoint));

% MAPEIBM1=sum(abs(real(ACOgBest.Error(:,1))))/(NumberOfTrainPoint);
% MAPEIBM2=sum(abs(real(testError(:,1))))/(NumberOfTestPoint);
% MAPEAPPLE1=sum(abs(imag(ACOgBest.Error(:,1))))/(NumberOfTrainPoint);
% MAPEAPPLE2=sum(abs(imag(testError(:,1))))/(NumberOfTestPoint);
% 
% MAPEDELL1=sum(abs(real(ACOgBest.Error(:,2))))/(NumberOfTrainPoint);
% MAPEDELL2=sum(abs(real(testError(:,2))))/(NumberOfTestPoint);
% MAPEMicroSoft1=sum(abs(imag(ACOgBest.Error(:,2))))/(NumberOfTrainPoint);
% MAPEMicroSoft2=sum(abs(imag(testError(:,2))))/(NumberOfTestPoint);
        
    time=toc;
    trial=['ACOResult_' ExperimentationName '_trial' int2str(trial)];
    save(trial);
    
 end