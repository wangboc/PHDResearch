function [Channel_1ToChannel_2, Channel_2ToChannel_1,coef_1, coef_2, ...
    GCSpectraChannel_1To_Channel_2, GCSpectraChannel_2To_Channel_1, ...
    PDCChannel_1ToChannel_2, PDCChannel_2ToChannel_1 ...
    RPCChannel_1ToChannel_2, RPCChannel_2ToChannel_1]  = ...
newFrequencyCausality(inChannel_1, inChannel_2, Fs, order)
% 参数说明：
% 输入参数说明：
% inChannel_1 = samples * trials; 一个通道的时序序列数据
% inChannel_2 = samples * trials; 另一个通道的时序序列数据
% order: 归模型的阶数
% Fs: 采样频率
% 
% 返回值参数说明：


% Channel_1ToChannel_2 = 1 * 1; 输入的第一个通道对第二个通道的因果关系值
% Channel_2ToChannel_1 = 1 * 1; 第二个通道对第一个通道的因果关系值
% coef_1 = (2*order) * 1: 用第一个通道的前order个数据序列加第二个通道的前order个数据序列
% 第一个通道的第order + 1个数据序列时前面（2*order）个序列点所对应的权重值（估计的平均值）
% coef_2 = (2*order) * 1： 用第二通道自身的的前order个数据序列，以及第一个数据通道的前order个
% 数据序列联合估计第二个通道的第order + 1个序列点时前2*order）个序列点所对应的权重值；

% 如果第一个类的trials不等于第二个通道的trials, 在调用本函数前请先处理

if size(inChannel_1, 2) ~= size(inChannel_2, 2)
    disp('输入参数不合法，请确认两个类的trials相等')
    return;
end

samples= size(inChannel_1, 1);
trials = size(inChannel_1, 2);

%% *********计算时域联合回归系数**********
%   coefficient_1 = (2*order) * trials 第一个通道的联合回归系数
%   coefficient_2 = (2*order) * trials 第二个通道的联合回归系数
% order = 3
for k = 1:trials  % 10
    [Ch_1_To_Ch_2, Ch_2_To_Ch_1, c1, c2] = ...
        newTimeCausality(inChannel_1(:, k), inChannel_2(:, k), order);
    coefficient_1(:, k) = c1;
    coefficient_2(:, k) = c2;
end

coef_1 = mean(coefficient_1, 2);
coef_2 = mean(coefficient_2, 2);

%% ---调整系数的排列顺序（血的教训，下次一定要测试好每一个结论）-----
coef_1_new(1:order,1) = coef_1(order:-1:1, 1);
coef_1_new(order+1:2*order,1) = coef_1(end:-1:order+1,1);

coef_2_new(1:order,1) = coef_2(order:-1:1, 1);
coef_2_new(order+1:2*order,1) = coef_2(end:-1:order+1,1);
%% 用bsmart里面计算系数的方法得出系数
% [A,Z] = compcoef([inChannel_1'; inChannel_2'],1,samples,order);
% jj = 1;
% for ii = 1 : 2 : order*2
%     coef_1_new(jj,1) = -A(1,ii);
%     coef_1_new(jj + 12) = -A(1,ii+1);
%     coef_2_new(jj,1) = -A(2,ii);
%     coef_2_new(jj + 12,1) = -A(2,ii+1);
%     jj = jj + 1;
% end
% coef_2_new = [coef_2_new(:,order+1 : 2*order);coef_2_new(:,1:order)]
%%
F = 0:0.05:30;
%% ***** 计算频域因果关系********
for k = 1:trials
    [delt_1_square_sum , delt_2_square_sum, co] = ...
    getDelt(inChannel_1(:, k), inChannel_2(:, k), coef_1, coef_2, order);
    [power1, fre_x] = pwelch(inChannel_1(:, k), 128, 120, F, Fs);
    [power2, fre_y] = pwelch(inChannel_2(:, k), 128, 120, F, Fs);
    
     power_x = power1';clear power1;
     power_y = power2';clear power2;
    % f0 = fre_x';
     %save fre f0
     fre = -2*pi*fre_x'/Fs;
     %fre = -2*pi*fre_x';
     cf11 = 0;
     cf12 = 0;
     cf21 = 0;
     cf22 = 0;
    
     for m = 1:order
         cf11 = cf11 + coef_1_new(m) * exp(i*m*fre);
         cf12 = cf12 + coef_1_new(order+m) * exp(i*m*fre);
         cf21 = cf21 + coef_2_new(m) * exp(i*fre*m);
         cf22 = cf22 + coef_2_new(order+m) * exp(i*fre*m);
     end
     
     a11f = abs(cf11);
     a12f = abs(cf12);
     a21f = abs(cf22);
     a22f = abs(cf21);
     %%  新因果关系
     %new_GCS_YX(k,:)=(a12f.*sqrt(power_y))./sqrt((a11f.^2.*power_x + a12f.^2.*power_y+delt_1_square_sum));
     %new_GCS_XY(k,:)=(a21f.*sqrt(power_x))./sqrt((a21f.^2.*power_x+a22f.^2.*power_y+delt_2_square_sum));
     new_GCS_YX(k,:)=(a12f.^2.*(power_y))./(a11f.^2.*(power_x) + a12f.^2.*(power_y)+delt_1_square_sum);
     new_GCS_XY(k,:)=(a21f.^2.*(power_x))./(a21f.^2.*(power_x)+a22f.^2.*(power_y)+delt_2_square_sum);
     
     %% 格兰杰频域因果关系
     GCa2 = 1-cf11;
     GCb2 = -cf12;
     GCd2 = -cf21; %注：由于系数排序的问题...
                   %因此对第二个通道来说cf21是第二个通道自己对自己的影响
     GCc2 = 1-cf22;%交换 GCd2 和  GCc2的值  下一个版本记得修复这个bug
     
     GCdetA = GCa2.*GCc2 - GCb2.*GCd2;
     Hxy = cf12./GCdetA;
     Hyx = cf22./GCdetA;    % 跟胡老师代码不同之处  还是因为系数排序问题
     Hxx = GCd2./GCdetA;
     Hyy = GCa2./GCdetA;
     
     S_xx1 = co(1,1)*abs((Hxx+(co(1,2)/co(1,1))*Hxy)).^2;
     S_xx2 = (co(2,2)-co(1,2)*co(1,2)/co(1,1))*abs(Hxy).^2;
     S_yy1 = co(2,2)*abs((Hyy+(co(2,1)/co(2,2))*Hyx)).^2;
     S_yy2 = (co(1,1)-co(1,2)*co(1,2)/co(2,2))*abs(Hyx).^2;
    
     GCS_YX(k,:)=log((S_xx1+S_xx2)./S_xx1);
     GCS_XY(k,:)=log((S_yy1+S_yy2)./S_yy1);
     PDC_YX(k,:)=abs(GCb2)./sqrt(abs(GCb2).^2+abs(GCd2).^2);
     PDC_XY(k,:)=abs(GCc2)./sqrt(abs(GCc2).^2+abs(GCa2).^2);
     RPC_YX(k,:)=abs(Hxy).^2*co(2,2)./(abs(Hxx).^2*co(1,1)+abs(Hxy).^2*co(2,2));
     RPC_XY(k,:)=abs(Hyx).^2*co(1,1)./(abs(Hyy).^2*co(2,2)+abs(Hyx).^2*co(1,1));
     
end

Channel_2ToChannel_1 = mean(new_GCS_YX, 1);
Channel_1ToChannel_2 = mean(new_GCS_XY, 1);

% granger
     GCSpectraChannel_2To_Channel_1 = mean(GCS_YX, 1);
     GCSpectraChannel_1To_Channel_2 = mean(GCS_XY, 1);
     PDCChannel_2ToChannel_1 = mean(PDC_YX, 1);
     PDCChannel_1ToChannel_2 = mean(PDC_XY, 1);
     RPCChannel_2ToChannel_1 = mean(RPC_YX, 1);
     RPCChannel_1ToChannel_2 = mean(RPC_XY, 1);


end

