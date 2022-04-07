c = cell(1, 100);
d = pi;
save('a.mat', 'c', 'd');
% clear
% cc = load('a.mat', 'c');   % ==> cc.c : {1x100 cell}