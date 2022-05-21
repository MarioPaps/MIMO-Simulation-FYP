%example
A = rand(5);
B = magic(10);
save example.mat A B -v7.3;
clear A B

exampleObject = matfile('example.mat');
exampleObject.Properties.Writable = true;
[nrowsB,ncolsB] = size(exampleObject,'B');
for row = 1:nrowsB
  exampleObject.B(row,:) = row * exampleObject.B(row,:);
end
%% 
exampleObject.test=  exampleObject.B *  exampleObject.B';
%% 
exampleObject = matfile('x50.mat');
exampleObject.Properties.Writable = true;
exampleObject.newx= exampleObject.x';
%% 

exampleObject.covmat= exampleObject.x * exampleObject.newx;