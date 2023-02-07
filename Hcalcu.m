function [H,Hz,Wy]=Hcalcu(csj,deltsita,fai,e,E,c,P,ps,Dc,Dh,Rc1,Rc2,n,m,WY0,lh)
% 计算各点处水膜厚度
dbstop if error
H=zeros(n+1,m+1);PP=zeros(n+1,m+1);Hz=zeros(n+1,m+1);
WY0c=WY0;
Wy=WyCalcu(P,m,n,E,ps,Dc,WY0c,lh);
% Wy=zeros(n+1,m+1);
Dh=[Dh;Dh(1)]/c;
for i=1:n+1 %周向
    jd=csj+(i-1)*deltsita-fai;
    for j=1:m+1 %轴向
        PP(i,j)=0;
        H(i,j)=1+e*(cos(jd))+Wy(i,j)/c+PP(i,j)-Rc1(i,j)/c-Rc2(i,j)/c+Dh(i);%初始弹性变形量
%         if H(i,j)<=5*10^(-6)
%             H(i,j)=0;
% %             disp “外载荷超出最大承载能力，建议增大转速”
%         end
        if H(i,j)<=5*10^(-6)
            H(i,j)=0;
            Hz(i,j)=1;
        end
    end
end
havera=c*mean(H(:));
if isnan(havera)
    disp 膜厚为NaN，出错
end
end