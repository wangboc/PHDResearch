function [GCT1To2, GCT2To1,coef_1,coef_2,b,r,delt_1,delt_2] = GrangerCausalityTime(inChannel_1, inChannel_2, order)
 %% 
 % 输入参数说明：
% inChannel_1 = samples * 1; 一个通道的时序序列数据
% inChannel_2 = samples * 1; 另一个通道的时序序列数据
% order: 归模型的阶数
% 
% 返回值参数说明：
% Channel_1ToChannel_2 = 1 * 1; 输入的第一个通道对第二个通道的因果关系值
% Channel_2ToChannel_1 = 1 * 1; 第二个通道对第一个通道的因果关系值
% coef_1 = (2*order) * 1: 用第一个通道的前order个数据序列加第二个通道的前order个数据序列估计
% 第一个通道的第order + 1个数据序列时前面（2*order）个序列点所对应的权重值
% coef_2 = (2*order) * 1： 用第二通道自身的的前order个数据序列，以及第一个数据通道的前order个
% 数据序列联合估计第二个通道的第order + 1个序列点时前2*order）个序列点所对应的权重值；

samples = size(inChannel_1, 1);
inChannel_1 = inChannel_1';
inChannel_2 = inChannel_2';
for i = (order+1):samples
    M1(i - order, :) = [inChannel_1(i-order:i-1), inChannel_2(i-order:i-1)];%M1  计算coef_1时用到的矩阵
    
    M2(i - order, :) = [inChannel_2(i-order:i-1), inChannel_1(i-order:i-1)];
end
% %***调试
% M1
% M7

y1 = inChannel_1(order+1:end)';
y2 = inChannel_2(order+1:end)';
 %[b,bint,r] = regress(y1,M1);
%   用最小二乘法求联合回归的系数
coef_1 = inv(M1'*M1) * (M1'*y1);
coef_2 = inv(M2'*M2) * (M2'*y2);
delt_1 = M1*coef_1 - inChannel_1(order+1:end)';
delt_2 = M2*coef_2 - inChannel_2(order+1:end)';

%delt_1_square_sum = delt_1' * delt_1
%delt_2_square_sum = delt_2' * delt_2


averagedelt_1 = mean(delt_1);
averagedelt_2 = mean(delt_2);
delt_1_square_sum = (delt_1 - averagedelt_1 )' * (delt_1 - averagedelt_1 ) / (samples - order - 1);%无偏估计方差
delt_2_square_sum = (delt_2 - averagedelt_2 )' * (delt_2 - averagedelt_2 ) / (samples - order - 1);
%% 自回归模型
for i = (order+1):samples
    M1_a(i - order, :) = inChannel_1(i-order:i-1);
    M2_a(i - order, :) = inChannel_2(i-order:i-1);
end
y1_a = inChannel_1(order+1:end)';
y2_a = inChannel_2(order+1:end)';
coef_1a = inv(M1_a'*M1_a) * (M1_a'*y1_a);
coef_2a = inv(M2_a'*M2_a) * (M2_a'*y2_a);
delt_1a = M1_a*coef_1a - inChannel_1(order+1:end)';
delt_2a = M2_a*coef_2a - inChannel_2(order+1:end)';
averagedelt_1a = mean(delt_1a);
averagedelt_2a = mean(delt_2a);
delt_1_square_suma = (delt_1a - averagedelt_1a )' * (delt_1a - averagedelt_1a ) / (samples - order - 1);
delt_2_square_suma= (delt_2a - averagedelt_2a )' * (delt_2a - averagedelt_2a ) / (samples - order - 1);
%%
GCT2To1 = log(var(delt_1a) / var(delt_1));
GCT1To2 = log(var(delt_2a) / var(delt_2));

