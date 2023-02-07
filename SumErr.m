function [sumt,sum1,sum2]=SumErr(PK,P,n,m)
%% Îó²î¿ØÖÆ
sum1=0;
sum2=0;
for s=1:n+1
    for t=1:m+1
        sum1=sum1+abs((PK(s,t)-P(s,t)));
        sum2=sum2+abs(P(s,t));
    end
end
sumt=sum1/sum2;
end