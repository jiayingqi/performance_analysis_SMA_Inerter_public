function [lamda, Phi, r]=complex_modes(M,C,K,E)
%****************************************************************
%----             Author(s): Wang Chao, Jia Yingqi           ----
%----             Affiliation: Tongji University             ----
%----             E-mail: jiayingqi@tongji.edu.cn            ----
%----             Date: 10/22/2020                           ----
%****************************************************************

[~,N]= size(M); % 自由度数
A=[zeros(N,N),M;M,C];
B=[-M,zeros(N,N);zeros(N,N),K];

F=[zeros(N,1);-M*E];
[Phi,lamda]=eig(-B, A);
% [Phi上半部分乘以过lamda，且以共轭形式出现，应调整顺序

lamda=diag(lamda);

% 归一化
Sum=conj(Phi).*Phi;
Phi=Phi./sqrt(sum(Sum,1));
% python程序中的eig特征向量是经过归一化处理的，对于matlab，应将每一列上的
% 元素均除以该列所有元素的模方的和的平方根（对于实数为均方根），对于复数而言，
% 模方等于复数本身乘以其共轭复数。

A1=Phi.'*A*Phi; %A1是非对角矩阵，此处复数矩阵的非共轭转置应用“.'”表示，
% “'”为求共轭（转置）矩阵
AD=diag(A1); %输出向量
r=Phi.'*F./AD; %等效和荷载向量