% H = GUI;
close all;

filename = '0185.jpg';
%filename = 'f0105.jpg';
I = imread(filename);   
 
imshow(I);
if size(I,3) == 3;
    I = rgb2gray(I);
end

outImg = imgaussfilt(I, 1);
%figure;
%imshowpair(I, outImg, 'montage');

% slider_val = 2;
% SigmaS = floor(9*slider_val+1);
% slider_val = 4;
% SigmaR = floor(19*slider_val+1);
% 
% tol = 0.01;
% Img=I;
% 
% % make odd
% if (mod(SigmaS,2) == 0)
%   w  = SigmaS + 1;
% else
%   w  = SigmaS;
% end
% [outImg, param] =  shiftableBF(double(Img), SigmaS, SigmaR, w, tol);
% 
% imshow(outImg);

strel_size1 = 9;
L = marker_watershed(outImg, strel_size1);
Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
figure;
imshow(I,[]);
hold on;
himage = imshow(Lrgb);
set(himage, 'AlphaData', 0.3);
hold off;
% imshow(Lrgb);
