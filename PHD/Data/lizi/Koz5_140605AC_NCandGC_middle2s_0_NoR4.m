clear all
load 140605/140605_Go_trials_data.mat
load 140605/140605_NoGo_trials_data.mat

Goeegnew=Goeeg(:,[1:28,30:48],509:2543); %去掉第29个trial (因为其全为0）
 NoGoeegnew=NoGoeeg(:,[1:5,7,9:48],509:2543);%去掉第6,8个trial (因为其全为0）
 clear Goeeg NoGoeeg;
 Goeeg=Goeegnew;
 NoGoeeg=NoGoeegnew;
  clear Goeegnew NoGoeegnew;
  
a=0;
order=10;
result1=[0,0];
result2=[0,0,0,0];
result1GC=[0,0];
result2GC=[0,0,0,0];
resultallNC=[0,0,0,0,0,0];
resultallGC=[0,0,0,0,0,0];
trialsGo=size(Goeeg,2);
trialsNoGo=size(NoGoeeg,2);
samples=size(Goeeg,3);
channels=size(Goeeg,1);
for i=1:channels
    for j=1:trialsGo
       
       GoEEG_R4mean=double(mean(Goeeg(i,j,:)));
       for h=1:samples
         Goeeg(i,j,h)=Goeeg(i,j,h)-GoEEG_R4mean;
       end
    end
end
for i=1:channels
    for j=1:trialsNoGo
       
       NoGoEEG_R4mean=double(mean(NoGoeeg(i,j,:)));
       for h=1:samples
         NoGoeeg(i,j,h)=NoGoeeg(i,j,h)-NoGoEEG_R4mean;
       end
    end
end
GoEEG_mean(1,:,:)=double(mean(Goeeg,1));
NoGoEEG_mean(1,:,:)=double(mean(NoGoeeg,1));
for h=1:1:channels
GoEEG_ave(h,:,:)=double(Goeeg(h,:,:))-a*GoEEG_mean;
NoGoEEG_ave(h,:,:)=double(NoGoeeg(h,:,:))-a*NoGoEEG_mean;
end;

for l=1:1:(channels-1)
   % for i=1:1:channels-1
    %for j=(i+1):1:channels
     for j=l+1:1:channels
         x(:,:)=GoEEG_ave(l,:,:);
         y(:,:)=GoEEG_ave(j,:,:);  
         x1(:,:)=NoGoEEG_ave(l,:,:);
         y1(:,:)=NoGoEEG_ave(j,:,:);  
       for k=1:1:trialsGo
           dataC3Go=x(k,:);dataC4Go=y(k,:);
           [chanel1(k),chanel2(k)] = newTimeCausality(dataC3Go',dataC4Go',order);
           [chanel1GC(k),chanel2GC(k)] = GrangerCausalityTime(dataC3Go',dataC4Go',order);
       end
       for k=1:1:trialsNoGo
           dataC3NoGo=x1(k,:);dataC4NoGo=y1(k,:);
           [chanel3(k),chanel4(k)] = newTimeCausality(dataC3NoGo',dataC4NoGo',order);
          [chanel3GC(k),chanel4GC(k)] = GrangerCausalityTime(dataC3NoGo',dataC4NoGo',order);
       end
           GoNewCausalityC3toC4=mean(chanel1);GoNewCausalityC4toC3=mean(chanel2);
          GoGCC3toC4=mean(chanel1GC);GoGCC4toC3=mean(chanel2GC);
           NoGoNewCausalityC3toC4=mean(chanel3);NoGoNewCausalityC4toC3=mean(chanel4);
          NoGoGCC3toC4=mean(chanel3GC);NoGoGCC4toC3=mean(chanel4GC);
           if 10000*(GoNewCausalityC3toC4-GoNewCausalityC4toC3)*(NoGoNewCausalityC3toC4-NoGoNewCausalityC4toC3)<0;
                               temp = [l,j],
        result1 = [result1;temp];
        result2 = [result2;[GoNewCausalityC3toC4,GoNewCausalityC4toC3,NoGoNewCausalityC3toC4,NoGoNewCausalityC4toC3]];
          size(result2,1)-1,
           end
            if 1000*(GoGCC3toC4-GoGCC4toC3)*(NoGoGCC3toC4-NoGoGCC4toC3)<0;
                  tempGC = [l,j],
        result1GC = [result1GC;tempGC];
        result2GC = [result2GC;[GoGCC3toC4,GoGCC4toC3,NoGoGCC3toC4,NoGoGCC4toC3]];
            end
          resultallNC=[resultallNC;[l,j,GoNewCausalityC3toC4,GoNewCausalityC4toC3,NoGoNewCausalityC3toC4,NoGoNewCausalityC4toC3]];
          resultallGC=[resultallGC;[l,j,GoGCC3toC4,GoGCC4toC3,NoGoGCC3toC4,NoGoGCC4toC3]];
     end
     [l,j],  
end
%result1,
result2,%size(result2),
%result1GC,
%result2GC,
%Causality_R5_R4.fig
%Causality_R5_L4.fig
save Koz5_140605resultNCandGC_AC_middle2s_0_NoR4 result1 result2 result1GC result2GC resultallNC resultallGC
%save 140715resultGC_AC result1GC result2GC
