function [M,C,K,E] = matrix_shear_building(m,c,k)
%****************************************************************
%----             Author(s): Wang Chao, Jia Yingqi           ----
%----             Affiliation: Tongji University             ----
%----             E-mail: jiayingqi@tongji.edu.cn            ----
%----             Date: 10/22/2020                           ----
%****************************************************************

zeta=0.05;
omg = sqrt(k./m);
N=length(m);
M=diag(m);

%% 阻尼矩阵
for i=1:N-1
    Cd(i,i)=c(i)+c(i+1);
    Cd(i,i+1)=-c(i+1);
    Cd(i+1,i)=Cd(i,i+1);
end
Cd(N,N)=c(N);

%% 刚度矩阵
for i=1:N-1
    K(i,i)=k(i)+k(i+1);
    K(i,i+1)=-k(i+1);
    K(i+1,i)=K(i,i+1);
end
K(N,N)=k(N);

%% 求频率和阻尼
% omega2 = eig(K, M);
% omega = sqrt(omega2);
% 
% w_i = omega(1);
% w_j = omega((floor(0.618*N))); % 模态截断
% 
% a0 = zeta*2.0*w_i*w_j/(w_i+w_j);
% a1 = zeta*2.0/(w_i+w_j);

% Cs = a0*M + a1*K % 瑞利阻尼
Cs = 0;
C = Cs + Cd; % 全部阻尼
E=ones(N,1);
