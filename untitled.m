% theta=[30;60]; %must be passed in as column vector
% 
% r=ones(3,5);
% S= r'*U;
%% loop
A = rand(3,7);
B = rand(3,7);
d = zeros(size(A,1)*size(B,1),1);
for i=1:size(A,2)
    d = d + kron(A(:,i),B(:,i));
end
%% no loop
dd = reshape(B*A',[],1);
%% 
ff = {[1 2; 3 4], [3 3; 1 1]}
fff=sum(cat(3,ff{:}),3)
%% 
A=ones(31,930);
out=squeeze(sum(reshape(A,31,30,[]),3));
%% 
A=ones(4,4);
dd=A*A';

B=A';
A1= A(1:height(A)/2,:);
A2=A(height(A)/2+1:end,:);
B1=B(:,1:width(B)/2);
B2=B(:,width(B)/2+1:end);
C=A1*B1+A2*B2;
