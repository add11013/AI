%Input is the data
%Output is IIM of data
function  output= CaculateIIM( data)

    for i=1:size(data,2)
     for j=1:size(data,2)
         
         %H(Y(X+)) ConditionEntropy: 1 is Positive,2 is Negative
         [h_yP,h_yOfxP]=ConditionEntropy(data,i,j,1);
         %H(Y(X-)) ConditionEntropy: 1 is Positive,2 is Negative
         [h_yN,h_yOfxN]=ConditionEntropy(data,i,j,2);
         
         xxx=data(:,i);
         pd=fitdist(xxx,'kernel');
         rat1=cdf(pd,0);
         rat2=1-rat1;
         if(i==j)         
            IIM(i,j)=0;
         else
            IIM(i,j)=(h_yP-h_yOfxP)*rat2+(h_yN-h_yOfxN)*rat1;
         end
     end
    end
 output=IIM;
    
end

