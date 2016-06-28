function [Channel_1ToChannel_2,Channel_2ToChannel_1,coef_1,coef_2]  = ...
newTimeCausality(inChannel_1, inChannel_2, order)

%      ,coef_1, coef_2, delt_1_square_sum, delt_2_square_sum
% ����˵����
% �������˵����
% inChannel_1 = samples * 1; һ��ͨ����ʱ����������
% inChannel_2 = samples * 1; ��һ��ͨ����ʱ����������
% order: ��ģ�͵Ľ���
% coef_1:ϵ����ֵ
% delt_1_square_sum:�в�����������ռ�ñ���
% 
% ����ֵ����˵����
% Channel_1ToChannel_2 = 1 * 1; ����ĵ�һ��ͨ���Եڶ���ͨ���������ϵֵ
% Channel_2ToChannel_1 = 1 * 1; �ڶ���ͨ���Ե�һ��ͨ���������ϵֵ
% coef_1 = (2*order) * 1: �õ�һ��ͨ����ǰorder���������мӵڶ���ͨ����ǰorder���������й���
% ��һ��ͨ���ĵ�order + 1����������ʱǰ�棨2*order�������е�����Ӧ��Ȩ��ֵ
% coef_2 = (2*order) * 1�� �õڶ�ͨ������ĵ�ǰorder���������У��Լ���һ������ͨ����ǰorder��
% �����������Ϲ��Ƶڶ���ͨ���ĵ�order + 1�����е�ʱǰ2*order�������е�����Ӧ��Ȩ��ֵ��

samples = size(inChannel_1, 1);
inChannel_1 = inChannel_1';
inChannel_2 = inChannel_2';
for i = (order+1):samples
    M1(i - order, :) = [inChannel_1(i-order:i-1), inChannel_2(i-order:i-1)];%M1  ����coef_1ʱ�õ��ľ���
    M2(i - order, :) = [inChannel_2(i-order:i-1), inChannel_1(i-order:i-1)];
end
% %***����
% M1
% M7
y1 = inChannel_1(order+1:end)';
y2 = inChannel_2(order+1:end)';

%   ����С���˷������ϻع��ϵ��
coef_1 = inv(M1'*M1) * (M1'*y1);
coef_2 = inv(M2'*M2) * (M2'*y2);
delt_1 = M1*coef_1 - inChannel_1(order+1:end)';
delt_2 = M2*coef_2 - inChannel_2(order+1:end)';

%delt_1_square_sum = delt_1' * delt_1
%delt_2_square_sum = delt_2' * delt_2

co = cov(delt_1, delt_2);  % ����ֵ��ʵ��ֵ֮������ƽ����
averagedelt_1 = mean(delt_1);
averagedelt_2 = mean(delt_2); 
%����ƫ�������ֱ����var
delt_1_square_sum = (delt_1 - averagedelt_1 )' * (delt_1 - averagedelt_1 ) / (samples - order - 1);
delt_2_square_sum = (delt_2 - averagedelt_2 )' * (delt_2 - averagedelt_2 ) / (samples - order - 1);
a = zeros(2);

%% ------����ڶ���ͨ���Ե�һ��ͨ����Ӱ��ֵ--------
M1_2 = M1(:, order+1:end);
M1_1 = M1(:, 1:order);
Channel_2InChannel_1Part = M1_2 * coef_1(order+1:end);  %�ڶ���ͨ���ڵ�һ��ͨ����ռ�ķ���
Channel_2InChannel_1PartSquareSum = ...
    Channel_2InChannel_1Part' * Channel_2InChannel_1Part; %�ڶ���ͨ���ڵ�һ��ͨ����ռ�ķ�����ƽ����
Channel_1InChannel_1Part = M1_1 * coef_1(1:order);  %��һ��ͨ���ڵ�һ��ͨ����ռ�ķ���

Channel_1InChannel_1PartSquareSum = Channel_1InChannel_1Part' * Channel_1InChannel_1Part;

%% ------�����һ����ͨ���Եڶ���ͨ����Ӱ��ֵ------
M2_1 = M2(:, order+1:end);
M2_2 = M2(:, 1:order);

Channel_1InChannel_2Part = M2_1 * coef_2(order+1:end);
Channel_1InChannel_2PartSquareSum = Channel_1InChannel_2Part' * Channel_1InChannel_2Part;
Channel_2InChannel_2Part = M2_2*coef_2(1:order); %������Լ���Ӱ ��
Channel_2InChannel_2PartSquareSum = Channel_2InChannel_2Part' * Channel_2InChannel_2Part;


% M2*coef_2
Channel_1ToChannel_2 = Channel_1InChannel_2PartSquareSum/(Channel_1InChannel_2PartSquareSum + Channel_2InChannel_2PartSquareSum + (samples)*var(delt_2));
Channel_2ToChannel_1 = Channel_2InChannel_1PartSquareSum/(Channel_1InChannel_1PartSquareSum + Channel_2InChannel_1PartSquareSum + (samples)*var(delt_1));
end