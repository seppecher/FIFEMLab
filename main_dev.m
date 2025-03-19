clear; close all;

I = double(imread('mediteranee.png'));
I = sum(I,3)/3;
I(I>0) = 1;

imagesc(I);
axis image
colormap gray

domain = Domain('image',I);

domain.plot('edgelabels');











