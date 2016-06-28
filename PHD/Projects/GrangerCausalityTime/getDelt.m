function [delt_1_square_sum , delt_2_square_sum, co] = ...
    getDelt(inChannel_1, inChannel_2, coef_1, coef_2, order)
% 参数说明：
% 输入参数说明：
% inChannel_1 = samples * 1; 一个通道的时序序列数据
% inChannel_2 = samples * 1; 另一个通道的时序序列数据
% order: 归模型的阶数
% coef_1: 第一个通道的回归模型平均系数
% coef_2: 第二个通道的回归模型平均系数
%
% 返回值参数说明：
% delt_1_square_sum = 1 * 1 用第一个联合回归模型系数估计时的误差平方的和
% delt_2_square_sum = 1 * 1 用coef_2估计第二个通道的数据时的误差平方的和

samples = size(inChannel_1, 1);
inChannel_1 = inChannel_1';
inChannel_2 = inChannel_2';
for i = (order+1):samples
    M1(i - order, :) = [inChannel_1(i-order:i-1), inChannel_2(i-order:i-1)];%M1  计算coef_1时用到的矩阵
    M2(i - order, :) = [inChannel_2(i-order:i-1), inChannel_1(i-order:i-1)];
end
y1 = inChannel_1(order+1:end)';
y2 = inChannel_2(order+1:end)';

delt_1 = M1*coef_1 - inChannel_1(order+1:end)';
delt_2 = M2*coef_2 - inChannel_2(order+1:end)';
%co = cov(delt_1, delt_2);
%delt_1_square_sum = delt_1' * delt_1;
%delt_2_square_sum = delt_2' * delt_2;
%----
delt_1_mean = mean(delt_1);
delt_2_mean = mean(delt_2);
co = cov(delt_1, delt_2);
delt_1_square_sum = (delt_1-delt_1_mean)' * (delt_1-delt_1_mean);
delt_2_square_sum = (delt_2-delt_2_mean)' * (delt_2-delt_2_mean);
end