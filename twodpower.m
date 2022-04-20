% function[power]= twodpower(in)
%     [N,len]=size(in);
%     antennapowers=zeros(N,1);
%     for n=1:N
%         antennapowers(n)= (1/len)*in(n,:)* ctranspose(in(n,:));
%     end
%     power= mean(antennapowers);
% %     acc=0;
% %     for ind=1:len
% %         singlecolpow= (1/Nbar)*ctranspose(in(:,ind))*in(:,ind);
% %         acc=acc+singlecolpow;
% %     end
% % 
% %     power=acc;
%    % power= sum(sum(abs(in).^2))/ (N*len); this gives the same answer as
%    % the loop
% end

function[power]=twodpower(in)
    [rows,cols]=size(in);
    powers=zeros(cols,1);
    for col=1:cols
        powers(col)=(1/rows)* in(:,col)'*in(:,col);
    end
    power=mean(abs(powers));

end