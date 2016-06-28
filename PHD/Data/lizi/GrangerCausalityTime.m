function [GCT1To2, GCT2To1,coef_1,coef_2,b,r,delt_1,delt_2] = GrangerCausalityTime(inChannel_1, inChannel_2, order)
 %% 
 % �������˵����
% inChannel_1 = samples * 1; һ��ͨ����ʱ����������
% inChannel_2 = samples * 1; ��һ��ͨ����ʱ����������
% order: ��ģ�͵Ľ���
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
 %[b,bint,r] = regress(y1,M1);
%   ����С���˷������ϻع��ϵ��
coef_1 = inv(M1'*M1) * (M1'*y1);
coef_2 = inv(M2'*M2) * (M2'*y2);
delt_1 = M1*coef_1 - inChannel_1(order+1:end)';
delt_2 = M2*coef_2 - inChannel_2(order+1:end)';

%delt_1_square_sum = delt_1' * delt_1
%delt_2_square_sum = delt_2' * delt_2


averagedelt_1 = mean(delt_1);
averagedelt_2 = mean(delt_2);
delt_1_square_sum = (delt_1 - averagedelt_1 )' * (delt_1 - averagedelt_1 ) / (samples - order - 1);%��ƫ���Ʒ���
delt_2_square_sum = (delt_2 - averagedelt_2 )' * (delt_2 - averagedelt_2 ) / (samples - order - 1);
%% �Իع�ģ��
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

