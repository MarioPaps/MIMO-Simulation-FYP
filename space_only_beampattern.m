%function to obtain figure beampatterns with space-only manifold vector

function[space_gain]= space_only_beampattern(DOAest,r,Fc,Fjvec,lightvel,Nc,Nsc)
    K=length(DOAest);
    N=length(r);
    space_gain=zeros(K*Nc*Nsc,360);
    for k=1:K
        acc=zeros(N,Nsc);
        for j=1:Nsc
           acc(:,j)=computeManifoldRx(DOAest(k),0,r,Fc,Fjvec(j),lightvel);
        end
        w= acc(:,1);

        for delk=0:Nc*Nsc-1
            pos=(delk+1)+(k-1)*Nc*Nsc;
            for theta=1:360
               space_gain(pos,theta)=w'*computeManifoldRx(theta,0,r,Fc,Fjvec(j),lightvel);
            end

        end
        surf((1:360),(0:Nc*Nsc-1),abs(space_gain((k-1)*Nc*Nsc+1:k*Nc*Nsc,:)),'FaceAlpha',1,'EdgeAlpha',0.5);
        hold on; %the hold on must always be here-otherwise the plot is not going to be correct
        xlabel('DOA (degrees)'); ylabel('Delay (Ts s)'); zlabel('Array Gain');
        shading('interp');
        colormap('jet');
   end

end






