function [ ]=SMA_SDOF_ST()
%% 本函数从频域和时域角度分别计算SMA-Inerter的响应
%****************************************************************
%----             Author(s): Wang Chao, Jia Yingqi           ----
%----             Affiliation: Tongji University             ----
%----             E-mail: jiayingqi@tongji.edu.cn            ----
%----             Date: 10/22/2020                           ----
%****************************************************************

clear
clc

%% 输入结构信息
m=20e3; % 原结构质量，kg
Tp=0.54; % s
omega=2*pi/Tp; % 原结构频率，rad/s
k=m*omega^2; % 原结构刚度，N/m
ksi=0.02; % 原结构阻尼比
c=2*ksi*omega*m; %原结构阻尼系数，N·s/m

mu=0.08; % 惯质比
kpas=0.16; % SMA element与原结构刚度之比
kpad=0.08; % 线性支撑弹簧与原结构刚度之比

xd1=0.005; % m
xd2=0.02; % m
ks=kpas*k; % SMA的初始刚度 N/m
alpha=0.001; % SMA屈服后和屈服前的刚度比

%% 无控结构
[M,C,K,E]=matrix_shear_building(m, c, k);
[lamda, Phi, r]=complex_modes(M,C,K,E);
[Omega, Sx1, Sigma_X1, Sigma_XP1, Sxp1, Sigma_XA1, Sxa1]=stochastic_response(lamda, Phi, r);

%% SMA控制的结构（ST1）
ke=100;
ce=50;
gap_k=ke;
gap_c=ce;
error=10^(-6);
% 迭代求解ke,ce，直到满足精度
while rms(gap_k)>error & rms(gap_c)>error
    temp_k=ke;
    temp_c=ce;
    kd_=alpha.*ks+(1-alpha).*ks.*ke;
    cd_=(1-alpha).*ks.*ce;
    [lamda, Phi, r] = complex_modes(m,c+cd_,k+kd_,1);
    [Omega, Sx2, Sigma_X2, Sigma_XP2, Sxp2, Sigma_XA2, Sxa2] = stochastic_response(lamda, Phi, r);
    ke=(xd2+xd1)./sqrt(2*pi)./Sigma_X2.*exp(1).^(-xd1.^2./2./Sigma_X2.^2);
    ce=(xd2-xd1)./sqrt(2*pi)./Sigma_XP2.*(1-erf(xd1/sqrt(2)./Sigma_X2));
    gap_k=ke-temp_k;
    gap_c=ce-temp_c;
end

%% SDI控制的结构(惯容与SMA先并联，再与普通弹簧串联)（ST2）
kd=kpad*k;
ke=10;
ce=50;
gap_k=ke;
gap_c=ce;
error=10^(-6);
% 迭代求解ke,ce，直到满足精度
item=1;
while rms(gap_k)>error & rms(gap_c)>error & item<30
    temp_k=ke;
    temp_c=ce;
    kd_=alpha.*ks+(1-alpha).*ks.*ke;
    cd_=(1-alpha).*ks.*ce;
    M=[m 0;0 mu*m];
    C=[c 0;0 cd_];
    K=[k+kd -kd;-kd kd_+kd];
    E=[1 0]';
    [lamda, Phi, r] = complex_modes(M,C,K,E);
    [Omega, Sx3, Sigma_X3, Sigma_XP3, Sxp3, Sigma_XA3, Sxa3] = stochastic_response(lamda, Phi, r);
    Sigma_XSMA3=Sigma_X3(2);
    Sigma_XPSMA3=Sigma_XP3(2);
    ke=(xd2+xd1)./sqrt(2*pi)./Sigma_XSMA3.*exp(1).^(-xd1.^2./2./Sigma_XSMA3.^2);
    ce=(xd2-xd1)./sqrt(2*pi)./Sigma_XPSMA3.*(1-erf(xd1/sqrt(2)./Sigma_XSMA3));
    gap_k=ke-temp_k;
    gap_c=ce-temp_c;
    item=item+1;
end

%% 减震比
% 位移响应计算
Sigma_X1
Sigma_X2
Sigma_X3

% 能量计算
ed_SMA=(1-alpha)*ks*...
    (xd2-xd1)/sqrt(2*pi)*(1-erf(xd1/sqrt(2)/Sigma_XSMA3))*Sigma_XPSMA3
ed_ST=c*Sigma_XP3(1)^2
etotal=ed_SMA+ed_ST;

disp('位移减震比')
sprintf('ST1    :%.4f',Sigma_X2(1)/Sigma_X1)
sprintf('ST2    :%.4f',Sigma_X3(1)/Sigma_X1)
disp('加速度减震比')
sprintf('ST1    :%.4f',Sigma_XA2(1)/Sigma_XA1)
sprintf('ST2    :%.4f',Sigma_XA3(1)/Sigma_XA1)
disp('能量过滤比')
sprintf('ST     :%.4f',1-ed_SMA/etotal)

%% 读取激励
file1=textread('..\地震波数据\GM1-TH.txt', '' ,'headerlines',1);
file2=textread('..\地震波数据\GM2-TH.txt', '' , 'headerlines',1);
file3=textread('..\地震波数据\GM3-TH.txt', '' , 'headerlines',1);
file4=textread('..\地震波数据\Artificial_EQSignal-TH.txt', '' , 'headerlines',1);

% 将地震波数组存储在元胞（矩阵的矩阵）中
amp=9.8; % m/s^2
wave{1}=file1(:,2)*amp;
wave{2}=file2(:,2)*amp;
wave{3}=file3(:,2)*amp;
wave{4}=file4(:,2)*amp;
dt=0.005; % 时间间隔

for i=1:4
    wave{i}=wave{i}'; %行转列
    wave{i}=wave{i}(:); %归为一列
    wave{i}=wave{i}';
    n(i)=length(wave{i});
end

%% Newmark求解位移时程响应
for i=1:4
    [u1{i},du1{i},ddu1{i}] = Newmark_belta(wave{i},dt,n(i),m, c, k,1);
    [u2{i},du2{i},ddu2{i}] = Newmark_belta(wave{i},dt,n(i),m,c+cd_,k+kd_,1);
    [u3{i},du3{i},ddu3{i}] = Newmark_belta(wave{i},dt,n(i),M,C,K,[1,0]');
    t{i}=linspace(0.005,n(i)*0.005,n(i));
end

%% 绘图
disp('计算全部完成，开始绘图')
blue=[55 126 184]/256;
orange=[255 160 65]/256;
green=[44 160 44]/256;
pink=[255 91 78]/256;
purple=[184 135 195]/256;
gray=[164 160 155]/256;

close all
% PSDF
% figure(1)
% semilogy(Omega/2/pi,[Sx1(1,:);Sx2(1,:);Sx3(1,:)],'linewidth',4)
% set(title('PSDF of displacement'),'Fontname', 'Times New Roman','FontSize',15)
% set(xlabel('Frequency (Hz)'),'Fontname', 'Times New Roman','FontSize',15)
% set(ylabel('PSDF'),'Fontname', 'Times New Roman','FontSize',15)
% set(legend('Original structure','SMA damper','SMA-Inerter'),'Fontname', 'Times New Roman','FontSize',15)
% set(gca,'Fontname', 'Times New Roman','FontSize',15)
% set(gcf,'position',[200,200,800,500])
% axis([0,2,10e-7,10e-2])
% set(gca,'looseInset',[0 0 0 0])
% grid on
% print('.\论文插图\Disp PSDF','-djpeg','-r200');
%
%
% 位移时程曲线
for i=1:4
    figure('position',[100,100,600,400])
    plot(t{i},u1{i}(1,:)*1e2,'linewidth',2,'color',blue)
    hold on
    plot(t{i},u2{i}(1,:)*1e2,'linewidth',2,'color',orange)
    plot(t{i},u3{i}(1,:)*1e2,'linewidth',2,'color',green)
    hold off
    if i==1
        ylim([-5,5])
        set(gca,'ytick',-5:2:5)
        set(gca,'yticklabel',sprintf('%.1f\n',get(gca,'ytick')))
    elseif i==2
        ylim([-2.5,2.5])
        set(gca,'ytick',-2.5:1:2.5)
        set(gca,'yticklabel',sprintf('%.1f\n',get(gca,'ytick')))
        xlim([0,25])
    elseif i==3
        ylim([-15,15])
        xlim([0,20])
        set(gca,'ytick',-15:5:15)
        set(gca,'yticklabel',sprintf('%.0f\n',get(gca,'ytick')))
    else
        ylim([-20,20])
        set(gca,'ytick',-20:5:20)
        set(gca,'yticklabel',sprintf('%.0f\n',get(gca,'ytick')))
    end
    set(xlabel('Time \it t \rm (s)'),'Fontname', 'Times New Roman','FontSize',15)
%     set(ylabel('Displacement (cm)'),'Fontname', 'Times New Roman','FontSize',15)
    set(ylabel('$$  \rm Displacement \; \it u_p \; \rm (cm) $$ ','interpreter','latex'),'Fontname', 'Times New Roman','FontSize',15)
    set(legend('Original structure','With conventional SMA damper','With SDI'),'Fontname', 'Times New Roman','FontSize',15,'EdgeColor',gray,'linewidth',1.5)
    set(gca,'Fontname', 'Times New Roman','FontSize',15,'linewidth',2)
    set(gca,'looseInset',[0 0 0 0],'linewidth',2,'Fontname','Times New Roman','FontSize',15)
    grid
    set(gca,'GridLineStyle', ':','GridColor','k')
%     print(['.\论文插图\Disp time history for GM ',num2str(i)],'-djpeg','-r300');
end
%
% 加速度
% for i=1:4
%     figure('position',[100,100,600,400])
%     plot(t{i},ddu1{i}(1,:),'linewidth',2,'color',blue)
%     hold on
%     plot(t{i},ddu2{i}(1,:),'linewidth',2,'color',orange)
%     plot(t{i},ddu3{i}(1,:),'linewidth',2,'color',green)
%     hold off
%     if i==1
%         ylim([-7.5,7.5])
%         set(gca,'ytick',-7.5:2.5:7.5)
%         set(gca,'yticklabel',sprintf('%.1f\n',get(gca,'ytick')))
%     elseif i==2
%         ylim([-3,3])
%         xlim([0,25])
%         set(gca,'ytick',-3:1:3)
%         set(gca,'yticklabel',sprintf('%.1f\n',get(gca,'ytick')))
%     elseif i==3
%         ylim([-20,20])
%         xlim([0,20])
%         set(gca,'ytick',-20:5:20)
%         set(gca,'yticklabel',sprintf('%.0f\n',get(gca,'ytick')))
%     else
%         ylim([-28,28])
%         set(gca,'ytick',-28:7:28)
%         set(gca,'yticklabel',sprintf('%.0f\n',get(gca,'ytick')))
%     end
%     set(xlabel('Time \it t \rm (s)'),'Fontname', 'Times New Roman','FontSize',15)
%     %     set(ylabel('Acceleration (m/s^2)'),'Fontname', 'Times New Roman','FontSize',15)
%     set(ylabel('$$  \rm Acceleration \; \it \ddot u_p \; \rm (m/s^2) $$ ','interpreter','latex'),'Fontname', 'Times New Roman','FontSize',15)
%     set(legend('Original structure','With conventional SMA damper','With SDI'),'Fontname', 'Times New Roman','FontSize',15,'EdgeColor',gray,'linewidth',1.5)
%     set(gca,'Fontname', 'Times New Roman','FontSize',15,'linewidth',2)
%     set(gca,'looseInset',[0 0 0 0],'linewidth',2,'Fontname','Times New Roman','FontSize',15)
%     grid
%     set(gca,'GridLineStyle', ':','GridColor','k')
%     print(['.\论文插图\Acc time history for GM ',num2str(i)],'-djpeg','-r300');
% end