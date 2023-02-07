function [PT,HT,Fx,Fy,Fp,Fc,Fmc,fmc]=RenoldRa(D,c,D2,L,lh,nd0,xs,w,E1,v1,W,DZ,DH)
dbstop if error
format long
tic
%%_________input variables_________________
% D  �ᾱֱ��
% c  �뾶��϶ ��λΪm
% L  ��г���
% l  �ڳĺ�2
% nd ˮ��ճ�� Pa.s
% w  ��תת�� r/min--rad/s
% E1  �ڳĵ���ģ��
% v1  �𽺲��ɱ�
% E2 �ᵯ��ģ��
% v2 �Ჴ�ɱ�
% v0 �������ɱ�
% W  ������غ�
% DZ ˮ������
% DH ˮ�۽Ƕ�
%_________output variable________
% H ���ˮĤ��ȷֲ�
% P ���ˮĤ
%___________testing data_________________
clc;close all;
% ps=1;DH=pi/50;DZ=5;
% D=0.361; c=0.0005; L=0.27;l=0.0145;nd=0.00089;w=20*2*pi/60; E1=870*10^6; v1=0.45; E2=207*10^9;v2=0.3; M=11.5; g=9.816;
% Qy=0;Qx=2000; Ra=5*10^(-6); rendaX=10*Ra;rendaY=10*Ra;Yp=24*10^6;gama=1;Hard=60*10^6;

ps=8000;DH=0;DZ=1;
D=0.236; c=0.0005; L=0.94;lh=0.0145;nd0=0.0008994;w=100*2*pi/60; E1=785*10^6; v1=0.47; E2=210*10^9;v2=0.3; M=11.5; g=9.816;
Qy=0;Qx=8000; Ra1=0.5*10^(-6); Rendax1=10*Ra1;Renday1=10*Ra1;Yp=24*10^6;gama=1;Hard=5*10^7;
Ra2=0.01*10^(-6);Rendax2=10*Ra2;Renday2=10*Ra2;
%__________formal function_________________
R=D/2;E=2*((1-v1^2)/E1+(1-v2^2)/E2)^(-1);
VV=6.*nd0.*w.*R.^2./(c.^2); %�ٶ������ٲ���
m=90; n=120; %�ᡢ����ȷ�����
Wx0=Qx+M*g;Fx0=0;
% ��z��һ������ z=ZR ����������ԳƷֲ�z=(-L/2, L/2)
deltL=L/m; %����ȷּ��
deltsita=(2*pi-DZ*DH)/DZ/n; %����Ƕȵȷִ�С
nk=round(DH/deltsita);Dr=deltsita;Dz=deltL;
e=0.9;fai=0; %��ʼƫ������ƫλ��
X0=c*e*cos(fai);Y0=c*e*sin(fai); %��ʼλ��
Xv0=0;Yv0=0; %�ٶȳ�ֵ
pbd=1;Dt=1/(10*Wx0/M);Dt0=Dt;%ʱ�������Dt<=1/���غ�
aT=10;ERR=0.001;ERRF=0.01;%���ƽ��λ��
XT=X0;YT=Y0; XvT=Xv0; YvT=Yv0;XaT=[];YaT=[];%�洢λ������
Dc=Dccal(m,n,Dr,Dz);
Rc1=Cucao(Ra1,Rendax1,Renday1,n,m,D,L,DZ);
Rc2=Cucao(Ra2,Rendax2,Renday2,n,m,D,L,DZ);
S=[];kk=1;faT=1;xs=1;
JHcs=[0.01,0.01];ks=1;
dyt=1;dxt=1;
while (abs(aT)>ERR)
    PT=[];HT=[];FN=[];WY=[];PPJ=[];fai0=fai;faT0=faT;
    for dk=1:DZ
        csj=(dk-1)*(2*pi/DZ);
        ccs=csj:deltsita:csj+(n-1)*deltsita;
        Dh=DhCalcu(R,ccs,DZ,JHcs,ks);
        ERR3=5.0e-3; %���
        %���㸳��ֵ����Ҫ��
        PD=zeros(nk,m);HD=20*ones(nk,m);Fn=zeros(n+1,m+1);Ar=0;Hr=0;
        %         cccp=10^(-10)*ones(n-1,m-1);PK=padarray(cccp,[1 1]);
        PK=zeros(n+1,m+1);PJ=zeros(n+1,m+1);WY0=zeros(n+1,m+1);
        k=1; %��������
        while k>0
            P=PK; PQ=PJ;%�´ε�����ֵ
            [H,Hz,Wy]=Hcalcu(csj,deltsita,fai,e,E,c,P,ps,Dc,Dh,Rc1,Rc2,n,m,WY0,lh);
            havera=c*mean(H(:));
            if isnan(havera)
                disp Ĥ��ΪNaN������
                return
            end
            Hxs=all(Hz);
            if Hxs==1
                for i=1:n+1
                    for j=1:m+1
                        if Hz(i,j)==1
                            P(i,j)=0;
                        end
                    end
                end
                PJ=Bianjie(Dc,ps,H,P,n,m,E,Hz);
                sumt2=abs((PQ-PJ)/PQ);
                if sumt2<=ERR3
                    break;
                end
            end
            PK2=P2calcu(H,Hz,deltsita,deltL,P,n,m,VV,c,e,w,R,csj,fai,ps);
            PK=PK+(PK2-PK)*xs;
            if k>10000
                break
            end
            [sumt,~,sum2]=SumErr(PK,P,n,m);
            if (sum2<10^(-5))&&(k>10)
                break;
            end
            if sumt<=ERR3
                break;
            end
            k=k+1;
            WY0=Wy;
        end
        PPK=PK(1:n,1:m);HK=H(1:n,1:m);
        PpJ=PJ(1:n,1:m);
        PT=[PT;PPK];HT=[HT;HK];WY=[WY;Wy];PPJ=[PPJ;PpJ];
        AR(dk)=Ar;HR(dk)=Hr;
    end
    PcxT=PT+PPJ;
    [Fx,Fy,~,Xa,Ya,Xv,Yv,X,Y,e,~,fai]=Fcalcu(deltsita,ps,PcxT,FN,DZ,R,L,c,Qx,Qy,M,g,Dt,X0,Y0,Xv0,Yv0,dyt,dxt,n,m);
    XT=[XT,X];YT=[YT,Y];XvT=[XvT,Xv];YvT=[YvT,Yv];XaT=[XaT,Xa];YaT=[YaT,Ya];
    X0=X;Y0=Y;%Xv0=Xv;Yv0=Yv;
    if abs(Fx)>Wx0/9
        Dt=0.3*Dt0;
    elseif abs(Fx)>Wx0/3
        Dt=0.01*Dt0;
    elseif abs(Fx)>Wx0*4/5
        Dt=0.001*Dt0;
    end
    if abs(Fx)<abs(Fx0)
        dxt=dxt/2;
    end
    Fx0=Fx;
    aT=abs((Wx0-Fx)/Wx0);
    pbd=pbd+1;
    if sum(sum(FN))~=0
        disp Ĥ��С�ڴֲڶȣ���������״̬
    end
    faT=100*abs((fai0-fai)/fai);
    S(pbd,:)=[Fx,Fy,Xa,Ya,X,Y,e,fai*180/pi,aT,faT];
end
% load DynamicLoad100rw2trttrpm100sin00025dt.mat
%Ħ��ϵ�����
toc
% [Fp,Fc,Fmc,fmc]=ucalcu(HT,PT,c,deltsita,deltL,R,w,Qx,Qy,M,g,nd0);
% ���Ĺ켣��ͼ
ZXplot(XT,YT,XvT,YvT,XaT,YaT);
%Ĥ���ͼ
PHplot(m,L,sitas,PT,HT)
end

%     aT=(Qy+M*g+Fy)/(M)
%     kk=abs(Fx/Fy)
%     if kk>0.2
%         fai=fai-(Fx/Fy)/(2*pi);
%     end
%     qq=atan(X/Y)



