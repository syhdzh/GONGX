function [Fm,Tm,Fn]=Propeller_weight(md,N,bmax,t02,t06,d,D,PD,W,rm,px,py)
%___________input variable__________________
% md ��Ҷ�ܶ�
% N ��Ҷ����
% bmax �����������
% t02 0.2R���������������
% t06 0.6R���������������
% d ���ֱ��
% D ������ֱ��
% PD ��������������������������յ�������,��λkw
% W ��������PD�µ�ת��,��λr/mi
% rm ������ҶƬ�����İ뾶
% px ����ƫ��x����
% py ����ƫ��y����

% Fn ��г��ܵľ����غ�
% Db ˮ������ھ�
% Lb ˮ����е����򳤶ȣ���
% au bu cu du Ħ��ϵ����ʽ�������
% vr ת�Ӻ�ˮ�����֮�����Ի����ٶ�
%-----------output variable-----------------
% Fm ��������������ƫ������ľ����غ�
% Tm ��������������ƫ�������Ť��
% Fnm ��ƽ���غɶ��������ѹ����Ӱ��
%______________testing data____________________
md=7600;N=3;bmax=0.132955656;
t02=0.155114932;t06=0.070909683;
d=0.886371041;D=4.3432181;
PD=850; W=100; rm=t02*0.25+t06*0.75; %����solidworks��������������ĵ�λ��
px=0.0001;py=0.0002;
au=0.0029;bu=4.621;cu=1.015;du=-0.3596;%������ˮ�������Ħ�����Լ����յ�����������ϵ���о�
%________________Weight caculation________________
ta=0.5*t02+t06;
da=(1-d/D)*D;
Gb=0.169*md*N*bmax*da*ta; %��Ҷ������,��λkg
Lk=d+0.1; %챳�����λm
K=1/13; %�����ϵ�׶��
d0=0.045+0.108*(PD/W)^(1/3)-K*Lk/2; %��챳������봦�ᾶ,m
Gm=(0.88-0.6*d0/d)*Lk*md*d^2; %�������,kg
Imp=9.48*md*N*bmax*ta*D^3; %������d/DС��0.18���������λkg`m`s2
%_______________weight center caculation____________________
len=360; %���ֵ�
a0=0:360/N:360*(N-1)/N; %����ҶƬ��ʼ�Ƕ�
x=zeros(len,1);y=zeros(len,1);rd=zeros(len,1); a=zeros(len,N);ar=zeros(len,1);
xa=zeros(len,N);ya=zeros(len,N);
for i=1:len
    for j=1:N
        a(i,j)=mod((a0(j)+(i-1)*360/len),360); %ҶƬ���ĽǶ�
        xa(i,j)=rm*cosd(a(i,j));
        ya(i,j)=rm*sind(a(i,j));
    end
    x(i)=sum(xa(i,:))/N; %����������x����
    y(i)=sum(ya(i,:))/N; %����������y����
    if abs(x(i))<1e-16
        x(i)=0;
    end
    if abs(y(i))<1e-16
        y(i)=0;
    end
    rd(i)=sqrt((x(i)-px)^2+(y(i)-py)^2); %����������ƫ����ת���ľ���
    if y(i)>=0
        ar(i)=atan(y(i)/x(i));
        if isnan(ar(i))
            ar(i)=pi/2;
        end
    else
        ar(i)=atan(y(i)/x(i))+pi;
        if isnan(ar(i))
            ar(i)=pi*3/2;
        end
    end
end
% scatter (x,y)
%______________________Fm calculation________________________
Fm=Gb+Gm;
Tm=Gb*rd+Gm*sqrt(px^2+py^2);
Fc=Gb*rd*(W*pi/30)^2;
alfaf=ones(len,1);
% alfaf��ת�����ذ�װ����ˮ���������֧��Ķ�̬�غɴ���Ч��
jd=(0:360/len:360*(len-1)/len)';
Fnm=alfaf.*Fc.*sind(jd);%��ƽ���غɶ��������ѹ����Ӱ��
plot(jd,Fn)
%��ƽ���غ������Ħ��ϵ���仯
P=Fn/(Db*Lb); %ԭ�غ��������ѹ��
PM=Fnm/(Db*Lb); %��ƽ���غ��������ѹ��
uf=au*exp(bu*(cu+abs(vr))^(-1))*(P+PM)^du;%���º��Ħ��ϵ��

end
