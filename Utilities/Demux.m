%inputs
%in: symbol stream for one user
%out: demultiplexed symbol stream for one user
function[out]= Demux(in,n,Nsc)
    out=zeros(Nsc, floor(n/Nsc));
    for ind=1:Nsc
          pos_start=1+(ind-1)*floor(n/Nsc);
          pos_end= ind*floor(n/Nsc);
          out(ind,:)= in(pos_start:pos_end);
    end
end