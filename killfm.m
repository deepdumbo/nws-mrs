function [y_sgsd,y_qz,y_fit]=killfm(y_ori,n)
% Algorithm to reduce symmetric frequency-modulated artifacts with opposite signals
% by extended-QZ algorithm
% Copyright (C) 2007-2012 Jyh-Miin Lin
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 2 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License along
%     with this program; if not, write to the Free Software Foundation, Inc.,
%     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
        % dlt=0.1;
dlt=0.5;
y_mag=real(y_ori).*(1+dlt)+i.*imag(y_ori).*(1-dlt);
% y_mag=real(y_ori).*2;
% for tmp_n=1;
%     n=256;
D=hankel(y_ori(1:2*n),y_ori(2*n:4*n));
M=hankel(y_mag(1:2*n),y_mag(2*n:4*n));
A=zeros(n,n,6);
A(:,:,1)=D(1:n,1:n);
A(:,:,2)=D(2:n+1,1:n);
A(:,:,3)=D(3:n+2,1:n);
A(:,:,4)=D(4:n+3,1:n);
A(:,:,5)=M(1:n,1:n);
A(:,:,6)=M(2:n+1,1:n);
A(:,:,7)=M(3:n+2,1:n);
A(:,:,8)=M(4:n+3,1:n);


A0=A;
[A2,Q,Z,V,W,A]=sgsd((A0),1,1);
% [A2,Q2,Z2,V2,W2,A]=sgsd(A0,1,1);
% A2=A;
ref_1=abs(diag(A(:,:,1)));

d_sgsd=(diag(A(:,:,2))./diag(A(:,:,1))+diag(A(:,:,3))./diag(A(:,:,2))+diag(A(:,:,4))./diag(A(:,:,3)))./3;
ki_sgsd_qz=((diag(A(:,:,5))./diag(A(:,:,1)))+(diag(A(:,:,6))./diag(A(:,:,2)))+ (diag(A(:,:,7))./diag(A(:,:,3)))+(diag(A(:,:,8))./diag(A(:,:,4))))./4;

ki_sgsd=((diag(A2(:,:,5))./diag(A2(:,:,1)))+(diag(A2(:,:,6))./diag(A2(:,:,2)))+ (diag(A2(:,:,7))./diag(A2(:,:,3)))+(diag(A2(:,:,8))./diag(A2(:,:,4))))./4;
% ki_sgsd=((diag(A(:,:,5))./diag(A(:,:,1)))+(diag(A(:,:,6))./diag(A(:,:,2)))+ (diag(A(:,:,7))./diag(A(:,:,3)))+(diag(A(:,:,8))./diag(A(:,:,4))))./4;
% thrld=(ki_sgsd<1.0  );

% [thrld d_sgsd ki_sgsd]
% [m0 m1 m2].'
for ind=1:n
    
    if (abs(ki_sgsd(ind))<3  && abs(d_sgsd(ind))<1)
    V(:,ind)=V(:,ind)./sqrt(V(:,ind).'*A0(:,:,1)*V(:,ind));
    aw1_sgsd(ind,1)=(A0(1,:,1)*V(:,ind)).^2;
para(ind,:)=[-real(log(d_sgsd(ind))) imag(log(d_sgsd(ind)))/2/pi abs(aw1_sgsd(ind,1)) imag(log(aw1_sgsd(ind,1)./abs(aw1_sgsd(ind,1))))];
para2(ind,:)=para(ind,:);
para3(ind,:)=para(ind,:);
tmp_ratio=((dlt+ki_sgsd(ind)-1)/dlt);
tmp_ratio_qz=((dlt+ki_sgsd_qz(ind)-1)/dlt);
if tmp_ratio>1
    tmp_ratio=1;
    tmp_ratio_qz=1;
end
para(ind,3)=tmp_ratio.*para(ind,3);
para3(ind,3)=tmp_ratio_qz.*para(ind,3);
    else
        aw1_sgsd(ind,1)=0;
%         d_sgsd(ind)=d_sgsd(ind)./abs(d_sgsd(ind));
para(ind,:)=zeros(1,4);
para2(ind,:)=para(ind,:);
para3(ind,:)=para(ind,:);

ki_sgsd(ind,:)=0;
ki_sgsd_qz(ind,:)=0;
    end
%     para(ind,3)=abs(aw1_qz(ind,1));
% para(ind,:)=[-real(d_sgsd(ind)) imag(d_sgsd(ind))/2/pi abs(aw1_sgsd(ind,1)) imag(log(aw1_sgsd(ind,1)./abs(aw1_sgsd(ind,1))))];
end
% para=para(para(:,3)>0,:);
% para2=para2(para(:,3)>0,:);
ki_sgsd=ki_sgsd(ki_sgsd>0,:);
ki_sgsd_qz=ki_sgsd_qz(ki_sgsd>0,:);
% para=itcmp(y_ori,-2)
% pole

y_sgsd=para2signal(para,y_ori);
y_fit=para2signal(para2,y_ori);
y_qz=para2signal(para3,y_ori);
return

function y2 = para2signal(para,y)
% Created by SRAM (Lin Jyh-miin Aug 17, 2006)
% Display result of exponential parameters of MRS
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
