function mainfunc (md,N,bmax,t02,t06,d,D,PD,W,rm)
%______________testing data____________________
md=7600;N=5;bmax=0.132955656;
t02=0.155114932;t06=0.070909683;
d=0.886371041;D=4.3432181;
PD=850; W=100; rm=t02*0.25+t06*0.75; %����solidworks��������������ĵ�λ��
px=0.0001;py=0.0002;
Qx=50;Qy=200; %������غ�
% ��о�ƽ��λ�ü���
[PT,HT,Fx,Fy,Fp,Fc,Fmc,fmc]=Reynold(D,c,D2,L,l,nd,xs,w,E,v,W,DZ,DH,Qx,Qy);
%���ƫ�������£�������صı仯
[Fm,Tm,Fn]=Propeller_weight(md,N,bmax,t02,t06,d,D,PD,W,rm,px,py);

