function [Fp,Fc,Fmc,fmc]=ucalcu(HT,PT,c,deltsita,deltL,R,w,Qx,Qy,M,g,nd,DZ,FN)
%% 求解摩擦力
HT0=HT;
HT0(isnan(HT0))=50*c;
[Fa,~]=gradient(PT); %求解压力的周向梯度
TFA=abs(Fa.*HT0);
Fp=sum(sum(TFA))*deltL*deltsita; %由液膜压力流产生的摩擦力
Fc=sum(sum(nd*(R+c)^2*w./HT0))*deltL*deltsita;%由液膜剪切流产生的摩擦力
Fmc=Fp+Fc; W=sqrt(Qx^2+(Qy+M*g)^2);
fmc=Fmc/W;
for i=1:n
    for j=1:m        
th(i,j)=nd*w*R/(HT(i,j)*c)+H(i,j)*c*PT(i+1,j)*ps/(deltsita*R);
tc(i,j)=uc*FN(i,j);


end
