%Marios Papadopoulos, EE4, 2022, Imperial College.
% 17/01/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takes the matrix of balanced gold codes and produces a matrix of balanced
% gold codes and a vector of the corresponding delays 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% in (Ncx(Nc+2))= matrix of gold codes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% bal_gold_codes = matrix of balanced gold codes (dimensions unkonwn
% a-priori)
%d elaypos= vector of delays for which balanced gold codes are obtained
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[bal_gold_codes,delaypos]=fFindBalancedGoldCodes(in)
    [rows,cols]=size(in);
    bal_gold_codes=[];
    delaypos=[];
    for j=1: cols
        num_ones= nnz(in(:,j));
        num_zeros= rows-num_ones;
        if(num_ones==num_zeros+1)
               bal_gold_codes=[bal_gold_codes in(:,j)];
               delaypos=[delaypos j];
        end
    end
end