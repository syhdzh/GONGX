function mainfunc (md,N,bmax,t02,t06,d,D,PD,W,rm)
%______________testing data____________________
md=7600;N=5;bmax=0.132955656;
t02=0.155114932;t06=0.070909683;
d=0.886371041;D=4.3432181;
PD=850; W=100; rm=t02*0.25+t06*0.75; %可用solidworks测算出螺旋桨质心的位置
px=0.0001;py=0.0002;
Qx=50;Qy=200; %轴承外载荷
% 轴承静平衡位置计算
[PT,HT,Fx,Fy,Fp,Fc,Fmc,fmc]=Reynold(D,c,D2,L,l,nd,xs,w,E,v,W,DZ,DH,Qx,Qy);
%轴承偏心作用下，轴承受载的变化
[Fm,Tm,Fn]=Propeller_weight(md,N,bmax,t02,t06,d,D,PD,W,rm,px,py);

