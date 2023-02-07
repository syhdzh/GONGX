function ZXplot(XT,YT,XvT,YvT,XaT,YaT)
%%画轴心轨迹图
scatter(XT,YT);xlabel('轴心轨迹x坐标(m)'); ylabel('轴心轨迹y坐标(m)')
saveas(gcf, '轴心轨迹位移坐标.jpg'); 
scatter(XvT,YvT);xlabel('轴心轨迹x速度(m/s)'); ylabel('轴心轨迹y速度(m/s)')
saveas(gcf, '轴心轨迹速度坐标.jpg');
scatter(XaT,YaT);xlabel('轴心轨迹x加速度(m2/s)'); ylabel('轴心轨迹y加速度(m2/s)')
saveas(gcf, '轴心轨迹加速度坐标.jpg');
TTT=(1:2001)/400;
plot(TTT,XT,TTT,YT);xlim([0 5]);xlabel('时间t(s)'); ylabel('轴心轨迹x,y位移(m)');legend('x位移','y位移')
saveas(gcf, '轴心x，y方向位移.jpg');
end