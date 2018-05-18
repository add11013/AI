clear
clc
OpenOriginalData=[];
OriginalData=[];

%% ----------------EX1----------------------

% data=xlsread('..\Data\EX1_TAIEX.csv');
% temp=data(~isnan(data(:,1)),1);
% OpenOriginalData=[OpenOriginalData temp];
% temp2=data(~isnan(data(:,4)),4);
% OriginalData=[OriginalData temp2];
% 
% data=xlsread('..\Data\EX1_DJIA.csv');
% temp=data(~isnan(data(:,1)),1);
% OpenOriginalData=[OpenOriginalData temp];
% temp2=data(~isnan(data(:,4)),4);
% OriginalData=[OriginalData temp2];
% 
% data=xlsread('..\Data\EX1_NASDAQ.csv');
% temp=data(~isnan(data(:,1)),1);
% OpenOriginalData=[OpenOriginalData temp];
% temp2=data(~isnan(data(:,4)),4);
% OriginalData=[OriginalData temp2];
% 
% data=xlsread('..\Data\EX1_SP500.csv');
% temp=data(~isnan(data(:,1)),1);
% OpenOriginalData=[OpenOriginalData temp];
% temp2=data(~isnan(data(:,4)),4);
% OriginalData=[OriginalData temp2];


%% ----------------EX2----------------------
data=xlsread('..\Data\EX2_TAIEX.csv');
temp=data(~isnan(data(:,1)),1);
OpenOriginalData=[OpenOriginalData temp];
temp2=data(~isnan(data(:,4)),4);
OriginalData=[OriginalData temp2];

data=xlsread('..\Data\EX2_HSI.csv');
temp=data(~isnan(data(:,1)),1);
OpenOriginalData=[OpenOriginalData temp];
temp2=data(~isnan(data(:,4)),4);
OriginalData=[OriginalData temp2];

% data=xlsread('..\Data\EX2_^N225.csv');
% temp=data(~isnan(data(:,1)),1);
% OpenOriginalData=[OpenOriginalData temp];
% temp2=data(~isnan(data(:,4)),4);
% OriginalData=[OriginalData temp2];
% 
% data=xlsread('..\Data\EX2_000001.SS.csv');
% temp=data(~isnan(data(:,1)),1);
% OpenOriginalData=[OpenOriginalData temp];
% temp2=data(~isnan(data(:,4)),4);
% OriginalData=[OriginalData temp2];

%% ----------------EX3----------------------
% data=xlsread('..\Data\EX3_APPLE.csv');
% temp=data(~isnan(data(:,1)),1);
% OpenOriginalData=[OpenOriginalData temp];
% temp2=data(~isnan(data(:,4)),4);
% OriginalData=[OriginalData temp2];
% 
% data=xlsread('..\Data\EX3_DELL.csv');
% temp=data(~isnan(data(:,1)),1);
% OpenOriginalData=[OpenOriginalData temp];
% temp2=data(~isnan(data(:,4)),4);
% OriginalData=[OriginalData temp2];
% 
% data=xlsread('..\Data\EX3_IBM.csv');
% temp=data(~isnan(data(:,1)),1);
% OpenOriginalData=[OpenOriginalData temp];
% temp2=data(~isnan(data(:,4)),4);
% OriginalData=[OriginalData temp2];
% 
% data=xlsread('..\Data\EX3_MS.csv');
% temp=data(~isnan(data(:,1)),1);
% OpenOriginalData=[OpenOriginalData temp];
% temp2=data(~isnan(data(:,4)),4);
% OriginalData=[OriginalData temp2];