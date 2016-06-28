function [delt_1_square_sum , delt_2_square_sum, co] = ...
    getDelt(inChannel_1, inChannel_2, coef_1, coef_2, order)
% ����˵����
% �������˵����
% inChannel_1 = samples * 1; һ��ͨ����ʱ����������
% inChannel_2 = samples * 1; ��һ��ͨ����ʱ����������
% order: ��ģ�͵Ľ���
% coef_1: ��һ��ͨ���Ļع�ģ��ƽ��ϵ��
% coef_2: �ڶ���ͨ���Ļع�ģ��ƽ��ϵ��
%
% ����ֵ����˵����
% delt_1_square_sum = 1 * 1 �õ�һ�����ϻع�ģ��ϵ������ʱ�����ƽ���ĺ�
% delt_2_square_sum = 1 * 1 ��coef_2���Ƶڶ���ͨ��������ʱ�����ƽ���ĺ�

samples = size(inChannel_1, 1);
inChannel_1 = inChannel_1';
inChannel_2 = inChannel_2';
for i = (order+1):samples
    M1(i - order, :) = [inChannel_1(i-order:i-1), inChannel_2(i-order:i-1)];%M1  ����coef_1ʱ�õ��ľ���
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