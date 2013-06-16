function [y,water] =nws2ws(x,M,N);
% matlab function using matrix-pencile method to generate WS data
% doing water suppression in postprocessing
% Copyright (C) 2007-2012 Jyh-Miin Lin except itcmp.m 
% itcmp.m was the work of its author Yung-Ya Lin
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
% input: x:NWS time domain fid. column vector
%           M: rank
%           N: signal length
%output:y:WS time domain fid, column vector
% usage e.g. y=nws2ws(x,7,1024);
% THIS function require itcmp.m and para2signal.m
if nargin==1
    M=15;
    N=round(length(x)/4);
%    N=512;
elseif nargin==2
     N=round(length(x)/4);
%    N=512;
else
end;
if N> length(x)
    N=length(x);
end
[para,tmpM]=itcmp(x(1:N),1);
para(1)=0;
para(3)=1;
x=x./para2signal(para,x);
% clear para tmpM;
[para,M]=itcmp((x(1:N)),M);
para=para(find(abs(para(:,2))<=0.05 | abs(abs(para(:,2))-0.5)<=0.05),:);
water=para2signal(para,x);
y=x-(water);


tmp_y=flipdim(y,1);
clear para;
[para,M]=itcmp(tmp_y(1:N),-2);
y=y-flipdim(para2signal(para,tmp_y),1);
return



function y2 = para2signal(para,y)
% Created by SRAM (Lin Jyh-miin Aug 17, 2006)
% Display result of exponential parameters of MRS
n=size(y);
[m,kk]=size(para);
    rec_s=complex(zeros(m,n(1)));
  for k=1:m
  rec_s(k,1:n(1))=para(k,3)*exp(1i*(para(k,4)))*exp((-para(k,1)+1i*2*pi*para(k,2))*((1:n(1))-1));
  end
  rec_s=permute(rec_s,[2 1]);
  y2=(sum(rec_s(:,:),2));
return
