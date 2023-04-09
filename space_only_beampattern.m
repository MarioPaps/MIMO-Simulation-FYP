%%MATLAB function to obtain the beampattern plot of the space-only MIMO beamformer
%DOAest: estimated DOAs of desired user
%r: Rx antenna array coordinates
%Fc: carrier frequency
%Fjvec: vector of subcarriers
%lightvel: velocity of light
%Nc: desired user's PN code length
%Nsc: number of subcarriers
%space_gain: beampattern matrix for all paths

function[space_gain]= space_only_beampattern(DOAest,r,Fc,Fjvec,lightvel,Nc,Nsc)
    K=length(DOAest);
    N=length(r);
    theta=(1:1:360);
    space_gain=zeros(K*Nc*Nsc,length(theta));
    for k=1:K
        acc=zeros(N,Nsc);
        accmanifolds=zeros(N,length(theta));
        for delk=0:Nc*Nsc-1
            pos=(delk+1)+(k-1)*Nc*Nsc;
            for j=1:Nsc
                    acc(:,j)=computeManifoldRx(DOAest(k),0,r,Fc,Fjvec(j),lightvel);
                    accmanifolds=accmanifolds+computeManifoldRx((theta)',0,r,Fc,Fjvec(j),lightvel);
            end
            %w=sum(acc,2); %using Nsc subcarriers multiplies the peak by Nsc
           
               w=acc(:,1);
               space_gain(pos,1:length(theta))=w'*computeManifoldRx((theta)',0,r,Fc,Fjvec(1),lightvel);
        end
        surf((theta),(0:Nc*Nsc-1),abs(space_gain((k-1)*Nc*Nsc+1:k*Nc*Nsc,:)),'FaceAlpha',1,'EdgeAlpha',0.5);
        hold on; %the hold on must always be here-otherwise the plot is not going to be correct
        xlabel('DOA (degrees)'); ylabel('Delay Ts (s)'); zlabel('Array Gain');
        xlim([0 380]);
        title('Massive MIMO Subspace Rx Beampattern');
        shading('interp');
        colormap('jet');
        ax = gca; 
        ax.FontSize = 11; 
    end

    %compute the 3 path weights, normalise, add them up and do the surf
    %plot
    %space_gain=zeros(Nc*Nsc,360);
    %     w=computeManifoldRx(DOAest,0,r,Fc,Fjvec(1),lightvel);
    %     w_norm= unitymag(w);
    %     w_tot= sum(w_norm,2,'omitnan');
    %     for delk=0:Nc*Nsc-1
    %             pos=(delk+1);
    %             space_gain(pos,1:360)=w_tot'*computeManifoldRx((1:360)',0,r,Fc,Fjvec(1),lightvel);
    %     end
    %     figure;
    %     surf((1:360),(0:Nc*Nsc-1),abs(space_gain),'FaceAlpha',1,'EdgeAlpha',0.5);
    %     xlabel('DOA (degrees)'); ylabel('Delay (Ts s)'); zlabel('Array Gain');
    %     shading('interp');
    %     colormap('jet');


end






