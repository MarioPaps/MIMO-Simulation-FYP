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
%% testbench
A=[1 2 ; 3 4; 5 6];
[out1,out2]= matrix_maxk(A,2);
%% convenient plotting
figure;
plot(ress,'-o');
ylim([-50 50]);
figure;
stem(ress,'o','LineStyle','none');
% x = 1:10;                                           % Create Data
% y = 2.5*rand(size(x));                              % Create Data
% stln = 0;
% figure
% plot(x, stln+zeros(size(x)), '-r')                 % Straight Line
% hold on
% plot(x, y, 'xk')                                    % Points
% plot([x; x], [zeros(size(y))+stln; y], '-g')        % Connecting Lines
% hold off
% ptlbls = compose('  %2.1f', y);
% text(x, y, ptlbls, 'HorizontalAlignment','left', 'VerticalAlignment','middle')
% axis([0  11    -0.1  2.5])

%% vectorisation hack
s=ones(9,3);
sb= 2*ones(3,1);
res= s .* sb'
%% betas randomly
R=0.5;
re= sqrt((R-0).*randn(1,10)+0).* cos(2*pi*rand(1,10));
re_abs= abs(re);
im= sqrt((R-0).*randn(1,10)+0).* sin(2*pi*rand(1,10));
im_abs= abs(im);

cond=1;
while(cond==1)
    indices= find(re_abs>R/2);
    if(isempty(indices))
        cond=2;
        break;
    end
    re(indices)=sqrt((R-0).*randn(1,length(indices))+0).* cos(2*pi*rand(1,length(indices)));
end


% beta= re+1i*im;
% norm(beta)
%% 
t=[3,4,5];
ee= find(t<2)
%% 
M=9;
N=9;
X = ones(M,N);
[M,N] = size(X);
Xcols = reshape(X,M,1,N);
Xrows = reshape(X,1,M,N);
result = bsxfun(@times,Xcols,Xrows)
%% more testing
x = findzero(2,3.5,0)
%% plot the antenna arrays
rx_bar=[2.56,2.37,1.81,0.98,0,-0.98,-1.81,-2.37,-2.56,-2.37,-1.81,-0.98,0,0.98,1.81,2.37]';
rx_bar=[rx_bar; rx_bar(1)];
ry_bar=[0,0.98,1.81,2.37,2.56,2.37,1.81,0.98,0,-0.98,-1.81,-2.37,-2.56,-2.37,-1.81,-0.98]';
ry_bar=[ry_bar; ry_bar(1)];
rz_bar=zeros(length(rx_bar),1);
figure;
plot(rx_bar,ry_bar,'.','MarkerSize',30,'DisplayName','Tx Antenna Array' );
rx=[1.43 0.90 -0.05 -0.98 -1.45 -1.24 -0.45 0.55 1.29]';
rx=[rx;rx(1)];
ry=[0.30 1.15 1.46 1.09 0.20 -0.77 -1.39 -1.36 -0.69]';
ry=[ry;ry(1)];
hold on;
plot(rx,ry,'.','MarkerSize',30,'DisplayName','Rx Antenna Array' );
grid on;
grid minor;
legend('show');
xlabel('x axis ({\lambda}/2)'); ylabel('y axis ({\lambda}/2)'); title('Tx and Rx Antenna Array Geometry');
ax = gca; 
ax.FontSize = 11; 
%% plot arrays
array1=[-2 0 0; -1 0 0; 0 0 0; 1 0 0; 2 0 0];
figure;
plot(array1(:,1),(zeros(1,5)),'-o');
% antenna array
increment= 360/5;
el_angles=[0,0+increment,0+2*increment,0+3*increment,0+4*increment];
R= 0.5/sind(increment/2);
rx= R*cosd(el_angles);
ry= R*sind(el_angles);
rz= R*zeros(size(rx));
array=[rx' ry' rz'];
figure;
plot(rx,ry,'.','MarkerSize',30,'DisplayName','5-element UCA' );
legend('show');

