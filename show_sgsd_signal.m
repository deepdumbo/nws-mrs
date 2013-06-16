function show_sgsd_signal
% demo function
m0=exp(-0.02+0.1i);
m1=exp(-0.02+0.2i);
m2=exp(-0.02+0.3i);
water0=exp(-0.05-0.4i);
t=(0:2047)';
db=40; 
itr=1;
% dlt=tmp_dlt/0.99-1;
dlt=0.5;

sigma=sqrt(0.1/(10^(db/10))); % noise amplitudes
rnd_noise=sigma*(randn(size(t))+i.*randn(size(t)));

% y_ori=m0.^t+m1.^t+m2.^t+i.*imag(water0.^t)+rnd_noise;
y0=m0.^t+m1.^t+m2.^t+1.*i.*imag(water0.^t);
y_ori=y0+rnd_noise;
% plot(real(fft(y_ori)))


%%%%%%%%%%%%%%%%%%%%%%%%core function
[y_sgsd y_qz y_fit]=killfm(double(y_ori),128);
%%%%%%%%%%%%%%%%%%%%%%%
% set(get(gcf,'CurrentAxes'),'FontName','Arial','FontSize',16);
% set(gcf,'DefaultLineLineWidth',1.5);
plot(real(fftshift(fft(y_ori))));hold on;
plot(real(fftshift(fft(y_fit))),'k');hold on;
plot(real(fftshift(fft(y_sgsd))),'r')
plot(real(fftshift(fft(y_qz))),'g')
legend('original signals','fitting of raw signals','FM removed by SGSD','FM removed by QZ')




% para
% size(y0)
% para(k,3)*exp(i*(para(k,4)))*exp((-para(k,1)+i*2*pi*para(k,2))*((1:n(1))-1));

% end
return

function y2 = para2signal(para,y)
% Created by SRAM (Lin Jyh-miin Aug 17, 2006)
% Display result of Matrix-pencil method (itcmp.m)
% please refer to itcmp.m
n=size(y);
[m,kk]=size(para);

if n(1)>=n(2)
  for k=1:m
  rec_s(k,1:n(1))=para(k,3)*exp(i*(para(k,4)))*exp((-para(k,1)+i*2*pi*para(k,2))*((1:n(1))-1));
  end
  rec_s=permute(rec_s,[2 1]);
  y2=(sum(rec_s(:,:),2));
elseif n(1)<=n(2)
  for k=1:m
  rec_s(k,1:n(2))=para(k,3)*exp(i*(para(k,4)))*exp((-para(k,1)+i*2*pi*para(k,2))*((1:n(2))-1));
  end
y2=(sum(rec_s(:,:),1));
end
return

