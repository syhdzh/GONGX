function [H,Hh,Ks]=Hcalcu2(deltsita,csj,fai,e,E,c,P,l,Ra,ps,Dc,Dh,n,m)
% 计算各点处水膜厚度
dbstop if error
H=zeros(n+1,m+1);Hh=zeros(n+1,m+1);PP=zeros(n+1,m+1);
Wy=zeros(n+1,m+1);
v0=0;nd=0.00089;w=100*2*pi/60;R=0.1805;
% Wy=WyCalcu(P,m,n,E,ps,Dc);
Dh=[Dh;Dh(1)]/c;
for i=1:n+1 %周向
    jd=(i-1)*deltsita+fai+csj;
    for j=1:m+1 %轴向
        PP(i,j)=6*l*((1-v0.^2)*nd*w*R.^2./(E*c.^3)).*P(i,j);
        H(i,j)=1+e*(cos(jd))+Wy(i,j)+Dh(i)+PP(i,j);%初始弹性变形量
        if H(i,j)<(Ra/c)
            Hh(i,j)=H(i,j);
        end
        if H(i,j)<0
            H(i,j)=0;
%             disp “外载荷超出最大承载能力，建议增大转速”
            continue
        end
    end
end
kp=sum(sum(Hh~=0));
Ks=kp/(m*n);
havera=l*mean(H(:));
if isnan(havera)
    disp 膜厚为NaN，出错
    return
end
end