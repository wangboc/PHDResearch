function [Channel_1ToChannel_2,Channel_2ToChannel_1,coef_1,coef_2]  = ...
newTimeCausality(inChannel_1, inChannel_2, order)

%      ,coef_1, coef_2, delt_1_square_sum, delt_2_square_sum
% 参数说明：
% 输入参数说明：
% inChannel_1 = samples * 1; 一个通道的时序序列数据
% inChannel_2 = samples * 1; 另一个通道的时序序列数据
% order: 归模型的阶数
% coef_1:系数的值
% delt_1_square_sum:残差项在整体中占得比重
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

%   用最小二乘法求联合回归的系数
coef_1 = inv(M1'*M1) * (M1'*y1);
coef_2 = inv(M2'*M2) * (M2'*y2);
delt_1 = M1*coef_1 - inChannel_1(order+1:end)';
delt_2 = M2*coef_2 - inChannel_2(order+1:end)';

%delt_1_square_sum = delt_1' * delt_1
%delt_2_square_sum = delt_2' * delt_2

co = cov(delt_1, delt_2);  % 估计值跟实际值之间误差的平方和
averagedelt_1 = mean(delt_1);
averagedelt_2 = mean(delt_2); 
%求无偏方差可以直接用var
delt_1_square_sum = (delt_1 - averagedelt_1 )' * (delt_1 - averagedelt_1 ) / (samples - order - 1);
delt_2_square_sum = (delt_2 - averagedelt_2 )' * (delt_2 - averagedelt_2 ) / (samples - order - 1);
a = zeros(2);

%% ------计算第二个通道对第一个通道的影响值--------
M1_2 = M1(:, order+1:end);
M1_1 = M1(:, 1:order);
Channel_2InChannel_1Part = M1_2 * coef_1(order+1:end);  %第二个通道在第一个通道中占的分量
Channel_2InChannel_1PartSquareSum = ...
    Channel_2InChannel_1Part' * Channel_2InChannel_1Part; %第二个通道在第一个通道中占的分量的平方和
Channel_1InChannel_1Part = M1_1 * coef_1(1:order);  %第一个通道在第一个通道中占的分量

Channel_1InChannel_1PartSquareSum = Channel_1InChannel_1Part' * Channel_1InChannel_1Part;

%% ------计算第一个个通道对第二个通道的影响值------
M2_1 = M2(:, order+1:end);
M2_2 = M2(:, 1:order);

Channel_1InChannel_2Part = M2_1 * coef_2(order+1:end);
Channel_1InChannel_2PartSquareSum = Channel_1InChannel_2Part' * Channel_1InChannel_2Part;
Channel_2InChannel_2Part = M2_2*coef_2(1:order); %自身对自己的影 响
Channel_2InChannel_2PartSquareSum = Channel_2InChannel_2Part' * Channel_2InChannel_2Part;


% M2*coef_2
Channel_1ToChannel_2 = Channel_1InChannel_2PartSquareSum/(Channel_1InChannel_2PartSquareSum + Channel_2InChannel_2PartSquareSum + (samples)*var(delt_2));
Channel_2ToChannel_1 = Channel_2InChannel_1PartSquareSum/(Channel_1InChannel_1PartSquareSum + Channel_2InChannel_1PartSquareSum + (samples)*var(delt_1));
end