function [A,Q,Z,V,W,A_qz]=sgsd(A,itr,ini)
% Simultaneously Decompose N-by-N-by-r (3D array) by extended-QZ algorithm
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
%% A: tensor, containing N-N-r high order array
%%    if r==2 the extended QZ is like conventional QZ algorithm
%%    such as [A(:,:,1), A(:,:,2)]=qz(A(:,:,1), A(:,:,2));
%% itr: times of iteration 
%% ini: started with QZ algorithm 
B=A;
r=size(A,3);
n=size(A,1);
Q=complex(eye(n));
Z=complex(eye(n));
V=Q;
W=Z;
% if r< 3
%     sprintf('the 3rd dimension of A must >= 3 !')
%    return 
% end

% if ini==1 %% initialization with qz
[A1,A2,Q,Z,V,W]=qz((A(:,:,1)),(A(:,:,2)));
A(:,:,1)=A1;
A(:,:,2)=A2;
% norm(Q)
  for running_r=1:r-2
    A(:,:,running_r+2)=Q*squeeze(A(:,:,running_r+2))*Z;
  end
% else
% end

A_qz=(A);

if (itr==0 && ini==1)
%     sprintf('Approximated results with QZ method')
    return
    
else  %% itr~=0, ini==1


 Q2=eye(n)+1;
 Z2=eye(n)+1;

 for pj=1:itr
% while ((norm(Q2-eye(n))-1)>eps )%||  norm(B(:,:,1)*V*diag(A(:,:,1)) - B(:,:,2)*V*diag(A(:,:,2)))>eps  )
     for j=1:n-1
           M=squeeze(A(j:n,j,:));
           [H,tmp1,tmp2]=svd(M);
%            s=diag(s);
%            s(2)/s(1)
%            Q2=complex(eye(n));
%            Q2(j:n,j:n)=H';
             q2=H';
         for running_r=1:r
            A(j:n,j:n,running_r)=q2*(A(j:n,j:n,running_r));
         end

      Q(j:n,j:n)=q2*Q(j:n,j:n);
     end  %%%% Q loop 

 
     for j=1:n-1
            N=(permute(A(n-j+1,1:n-j+1,:),[ 3,2,1 ]));
            Z2=complex(eye(n));
            [tmp1,tmp2,v]=svd(N);
%             G=flipdim(v,2);
%             G=v(:,end:-1:1);
%             Z2(1:n-j+1,1:n-j+1)=G;
             z2=v(:,end:-1:1);
          for running_r=1:r           
              A(1:n-j+1,1:n-j+1,running_r)=(A(1:n-j+1,1:n-j+1,running_r))*z2;
          end 
            Z(1:n-j+1,1:n-j+1)=Z(1:n-j+1,1:n-j+1)*z2;

      end  %%%% Z loop 
    
%       V=V*Z2;
%       W=(Q2'*W')';
%        W=W*Q2;
    
end %%%% itr loop 


% tmp_K=zeros(n,n);
 
%      for j=1:n %% invariant subspace
%           for running_r=1:r      
%               
%               tmp_K(:,:)=tmp_K-(-1)^(running_r).*B(:,:,running_r).*A(j,j,running_r);
%               
%           end 
% %          tmp_K=B(:,:,1).*A(j,j,1)- B(:,:,2).*A(j,j,2);
%             for pj=1:r
%                  
%               V(:,j)=tmp_K\B(:,:,r-pj+1)*V(:,j);
% %              V(:,j)=tmp_K\B(:,:,pj)*V(:,j);
%               V(:,j)=V(:,j)./norm(V(:,j));
%              end 
% 
%       end %% invariant subspace 

     
 
 
end  %% end if (itr==0, ini==1)
end
