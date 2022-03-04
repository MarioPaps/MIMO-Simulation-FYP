%Marios Papadopoulos, EE4, 2022, Imperial College.
% 17/01/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reads an image file with AxB pixels and produces a column vector of bits
% of length Q=AxBx3x8 where 3 represents the R, G and B matrices used to
% represent the image and 8 represents an 8 bit integer. If P>Q then
% the vector is padded at the bottom with zeros.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% filename (String) = The file name of the image
% P (Integer) = Number of bits to produce at the output - Should be greater
% than or equal to Q=AxBx3x8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% bitsOut (Px1 Integers) = P bits (1's and 0's) representing the image
% x (Integer) = Number of pixels in image in x dimension
% y (Integer) = Number of pixels in image in y dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bitsOut,x,y]=fImageSource(filename,P)
%     image= importdata(filename);
%     [x,y,~]= size(image);
%     %%convert image to binary
%     bitsOut=dec2bin(image);
%     binvec= logical(bitsOut - '0');
%     bitstream= reshape(binvec',P,1); %convert to column vector
%     bitstream= double(bitstream);
%     bitsOut= bitstream; %final bitstream output

    i = importdata(filename);
    images = imresize(i,1);
    [x,y,z] = size(images);
    bitsOut = dec2bin(images(:))';
    bitsOut = double(bitsOut(:))-48; 
    if(length(bitsOut)<P)
        bitsOut=[bitsOut;zeros(P-length(bitsOut),1)];
    end
 
end