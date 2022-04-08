%Marios Papadopoulos, EE4, 2022, Imperial College.
% 17/01/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the STAR manifold vector for particular delay and azimuth pair
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% delay (scalar)= delay value
% theta (scalar)= azimuth value
% J  = shifting matrix
% pncode (Wx1 Integer) = W values of 1's and -1's
% array (NXM complex) = antenna array geometry matrix  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% h (N*Next x 1 complex) = computed spatiotemporal manifold 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [h] = getSpatiotemporalManifold(delay,theta,J,pncode,array,Fc,Fj,lightvel)
    Nc=length(pncode);
    c= [pncode;zeros(length(J)-Nc,1)];
    S= computeManifoldRx(theta,0,array,Fc,Fj,lightvel);
    h= kron(S,((J)^delay)*c);
end