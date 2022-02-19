%Marios Papadopoulos, EE4, 2022, Imperial College.
% 17/01/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takes two M-Sequences of the same length and produces a gold sequence by
% adding a delay and performing modulo 2 addition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% mseq1 (Wx1 Integer) = First M-Sequence
% mseq2 (Wx1 Integer) = Second M-Sequence
% shift (Integer) = Number of chips to shift second M-Sequence to the right
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% GoldSeq (Wx1 Integer) = W bits of 1's and 0's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [GoldSeq]=fGoldSeq(mseq1,mseq2,shift)
     m_seq2_shifted= circshift(mseq2,shift);
     GoldSeq= bitxor(mseq1,m_seq2_shifted);
end