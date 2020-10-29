function [exponent,Amplitude,Alpha,time,st,epochs]=dfa_2018(signal,start,stop,num_segment,start_fit,stop_fit,prune,fg);
% Calculates DFA, time is in samples, step is a multiplying factor
% signal: is signal in samples
% start: length of first segment (e.g. 3 seconds, calculate in samples)
% stop: length of last segment e.g. 55 seconds (in samples)
% num_segm: number of segments (e.g. 20)
% start_fit as start
% stop_fit as stop
% prune: is two column vector with starting points and duration of good segments 
% Example
% % Example with vector V, from 3 to 50 seconds, sampling freq = 200, 30 segments 
% N=randn(20*60*200,1); [b,a]=butter(2,[10 12]/100);
% 
% Na=abs(hilbert(filtfilt(b,a,N)));
% prune=[1 length(Na)-1];
% [exponent,Amplitude,Alpha,time,st,epochs]=dfa_data_general_mod(Na,3*200,50*200,30,3*200,50*200,[1 length(Na)-1],1);
Prune_time_1=prune(:,1);
Prune_time_2=prune(:,2);


step=(log10(stop)-log10(start))/(num_segment-1); 
 KK=length([log10(start):step:log10(stop)+eps]);
      for KKn=1:KK
      eval(['Alpha' num2str(KKn) ' = []; ' ])
      end
 
Amplitude=[];

for prune_num=1:size(Prune_time_1,1)
 k=0;
 vect1=signal(Prune_time_1(prune_num):Prune_time_1(prune_num)+Prune_time_2(prune_num)); vect1=vect1(:); 
 
% Amplitude(prune_num)=mean(vect1.^2);
 Amplitude=[Amplitude; vect1]; 
 vect1=vect1-mean(vect1); % centre the signal
 vect1=cumsum(vect1); % cumulative sum
 step=(log10(stop)-log10(start))/(num_segment-1); 
   for nn=log10(start):step:log10(stop)+eps
   n=round(10^nn); 
   K=floor(length(vect1)/n);
   k=k+1;
   time(k)=n;
   vect1R=reshape(vect1(1:K*n),n,K);
   %%%%%%%%%%%%%%%% Additional part with flipping
   vect2=flipud(vect1);
   vect1R_2=reshape(vect2(1:K*n),n,K);
   vect1R=[vect1R vect1R_2];
   %%%%%%%%%%%%%
   vect1RD=detrend(vect1R);
   %vect1RD=reshape(vect1R,1,2*K*n);
   %alpha(prune_num,k)=(mean(vect1RD.^2));
   % eval(['alpha' num2str(k) ' = mean(vect1RD.^2);'])
   tmp1 = mean(vect1RD.^2);
   eval(['Alpha' num2str(k) ' = [tmp1 Alpha' num2str(k) ']; ' ])
   end
 end

%alpha=sqrt(mean(alpha,1));
for kn=1:k
eval(['Alpha(' num2str(kn) ') = sqrt(mean(Alpha' num2str(kn) '));'])
end

Amplitude=mean(Amplitude);
I1=find(time>=start_fit); I1=I1(1);
I2=find(time>=stop_fit); I2=I2(1);
X=[ones(I2-I1+1,1),log10(time(I1:I2))'];
Y=log10(Alpha(I1:I2))';
[slope,bint,r,rint,stats]=regress(Y,X);
exponent=slope(2,1);
conf=((bint(2,2))-(bint(2,1)))/2;	% compute +- 95% confidence intervals
st=stats(1,1);
%if stats(1,1)<0.95; 'Warning R^2 values are low', end
if fg==1
    figure
    plot(X(:,2),Y,'ro')
    lsline
    hold
    plot(log10(time),log10(Alpha),'.')
    axis([min(log10(time)) max(log10(time)) min(log10(Alpha)) max(log10(Alpha))])
    title(['DFA, 1 vector, mean Alpha exp=', num2str(slope(2,1),3), ' +-=', num2str(conf,3),', R^2 = ', num2str(stats(1,1),3)])
    xlabel('log10 n, (time window in samples)')
    ylabel('log10 F(n)')
end
size(time);
epochs=length(Prune_time_1);
%end






