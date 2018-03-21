%Input of this function is Original data
%Output of this function is feature index(selected) and the DataMatrix
function [output,DataMatrix]=FeatureSelection(OriginalData, NumberOfTrainPoint, ExperimentationName)
%% IIM
% get feature
%caculate Number Of Target
NumberOfTarget=size(OriginalData,2);
%cacurate diff. of all Feature and save as data(t).value
for t=1:NumberOfTarget
    DataTemp=OriginalData(:,t);
    ColOfTarget=31;%(���^)+1
    LengthOfData=length(DataTemp);
    for jj=1:ColOfTarget
        k=jj;
        for i=1:LengthOfData-ColOfTarget
            TMP(i,jj)=DataTemp(k+1)-DataTemp(k);
            k=k+1;
        end
    end
    data(t).value=TMP;
end

    
    %combine all Feature(no target)
    AllFeature=[ ];
    for i=1:NumberOfTarget
        TMP=data(i).value(:,1:ColOfTarget-1);
        AllFeature=[AllFeature TMP];
    end
    %caculate data matrix
    DataMatrix=AllFeature;
    for i=1:NumberOfTarget
        temp=data(i).value(:,ColOfTarget);
        DataMatrix=[DataMatrix temp];
    end
    
    %combine the target(AllFeature+1target)
    for i=1:NumberOfTarget
        llMData(i).value=[AllFeature data(i).value(:,ColOfTarget)];
    end
    
  %use function "CaculateIIM" to caculate all IIM for feature selection
%      for i=1:NumberOfTarget
%          IIM(i).value=CaculateIIM(llMData(i).value(1:NumberOfTrainPoint,:));
%          i
%      end
        load(['IIM_' ExperimentationName]);
    
%% caculate selection gain
%initiate each target's "selected feature pool"(SP)
for j=1:NumberOfTarget
    PreSP(j).value(1)=0;
end
%��X��Ugain�̤j���S�x
for j=1:NumberOfTarget
    %initiate index
    for i=1:length(IIM(j).value)-1
        index(i)=i;
    end
    count=1;
    %�n�Dlength(IIM(j).value)-1��(�C�ӯS�x���n�D)
    for ite=1:length(IIM(j).value)-1
        %��gain�ȵ�SPold(j).gain
        %�ѤU�S�Q�D�L�����n��gain�ȡA�ҥH��ite�^�X�n��length(IIM(j).value)-ite��
        for i=1:length(IIM(j).value)-ite
            Redundancy=0;
            i1=index(i);
            %�p�G���F��N�⤾�l��T�q�A�S�F��N������v�T��T�q��@gain
            if ~PreSP(j).value==0
                for ii=1:length(PreSP(j).value)
                    i2=PreSP(j).value(ii);
                    InformationFi1TOFi2=IIM(j).value(i1,i2);
                    InformationFi2TOFi1=IIM(j).value(i2,i1);
                    temp=(InformationFi1TOFi2+InformationFi2TOFi1);
                    Redundancy=Redundancy+temp;
                end
                Information=IIM(j).value(i1,size(IIM(j).value,2));
                gain=Information-(Redundancy/(2*length(PreSP(j).value)));
%                 gain=Information-Redundancy;
            else
                gain=IIM(j).value(i1,size(IIM(j).value,2));
            end
            SPold(j).gain(i)=gain;
        end
            
            %��̤j
            tag=1;
            i=1;
            while tag==1
                if  SPold(j).gain(i)==max(SPold(j).gain)
                    PreSP(j).gain(count)=SPold(j).gain(i);
                    SPold(j).gain(i)=[];
                    PreSP(j).value(count)=index(i);
                    index(i)=[];
                    count=count+1;
                    tag=0;
                end
                i=i+1;
            end
    end
end

%��PreSP�̭������Ȩ��X�ӵ�SP
for j=1:NumberOfTarget
    count=1;
    for i=1:length(PreSP(j).value)
        if PreSP(j).gain(i)>0
            SP(j).gain(count)=PreSP(j).gain(i);
            SP(j).value(count)=PreSP(j).value(i);
            count=count+1;
        end
    end
end

%% caculate Omega�AOmega���C��SP�̪��S�x���X(�Y�����ưO�@�ӴN�n)
    %�N�Ĥ@��SP�̪��S�x�����[�JOmega
    Omega=SP(1).value;
    LengthOfOmega=length(Omega);
    %�@�}�l�w�g�NSP(1)�����[�JOmega�ҥH�q2�}�l
    for i=2:NumberOfTarget
        for ii=1:length(SP(i).value)
            %any (x==a)�p�Gx�����@�өΦh��a�^��1�A�S�h�^��0
            %�YOmega���S��SP(i).value(ii)�A�NSP(i).value(ii)�[�JOmega
            if ~any(Omega==SP(i).value(ii))
                LengthOfOmega=LengthOfOmega+1;
                Omega(LengthOfOmega)=SP(i).value(ii);
            end
        end
    end
    %% caculate NOL�ANOL���S�x�ثe�X�{���ֿn����
    NOL(1:LengthOfOmega)=0;
    for ii=1:LengthOfOmega
        for i=1:NumberOfTarget
            %��Omega(ii)���S���bSP(i).value�̡A��NOL�N+1
            for iii=1:length(SP(i).value)
                if Omega(ii)==SP(i).value(iii)
                    NOL(ii)=NOL(ii)+1;
                end
            end
        end
    end

    %% caculate w�Aw is covering rate 
    for i=1:length(NOL)
        w(i)=NOL(i)/NumberOfTarget;
    end
    wMean=mean(w);

    %% caculate gsum�Agsum is sum of the feature's selection gain in each selected feature pool
    gSum(1:LengthOfOmega)=0;
    for i=1:length(Omega)
        for ii=1:NumberOfTarget
            for iii=1:length(SP(ii).value)
                if Omega(i)==SP(ii).value(iii)
                    gSum(i)=gSum(i)+SP(ii).gain(iii);
                end
            end
        end
    end
    gSumMean=mean(gSum);
    


    %% caculate p�Ap is contribution index which is covering rate by gsum (W*gsum)
    NumberOfFP=0;
    for i=1:LengthOfOmega
        p(i)=w(i)*gSum(i);
    end
    
    %% Sort p and Omega
    for i=1:length(p)
        for ii=1:length(p)
            if p(ii)<p(i)
                pTemp=p(i);
                p(i)=p(ii);
                p(ii)=pTemp;

                OmegaTemp=Omega(i);
                Omega(i)=Omega(ii);
                Omega(ii)=OmegaTemp;
            end
        end
    end
    
    FP.index=Omega;
    FP.p=p;

    %% set upper bound and the lower bound�Aavoid much or less features in final pool
    Upper=4;
    lower=2;

    %output is final selected feature index for multi-target
    if length(FP.index)>Upper
        output=FP.index(1:Upper);
    elseif length(FP.index)<lower
        FP.index=Omega(1:lower);
        output=FP.index;
    else
        output=FP.index;
    end
toc
end