function L = marker_watershed(img, strel_size1)
    I = uint8(img);

    se = strel('disk',strel_size1); %((((((((((((((param_to_be_tuned)))))))))))))))))))
    Ie = imerode(I, se);
    Iobr = imreconstruct(Ie, I);
    % figure, imshow(Iobr), title('Opening-by-reconstruction (Iobr)')
    Iobrd = imdilate(Iobr, se);
    Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
    Iobrcbr = imcomplement(Iobrcbr);
    % figure, imshow(Iobrcbr), title('Opening-closing by reconstruction (Iobrcbr)')

    % use Iobrcbr for marker extraction
    fgm = imregionalmax(Iobrcbr);
    % figure, imshow(fgm), title('Regional maxima of opening-closing by reconstruction(fgm)')
    I2 = I;
    I2(fgm) = 255;
    fgm4 = bwareaopen(fgm, strel_size1);
    I3 = I;
    I3(fgm4) = 255;
    % figure, imshow(I3), title('Modified regional maxima superimposed on original image (fgm4)')

    % Step 3: compute foreground markers
    bw = im2bw(Iobrcbr, graythresh(Iobrcbr));
    % figure, imshow(bw), title('Thresholded opening-closing by reconstruction (bw)')

    D = bwdist(bw);
    DL = watershed(D);
    bgm = DL == 0;
    % figure, imshow(bgm), title('Watershed ridge lines (bgm)')

    hy = fspecial('sobel');
    hx = hy';
    Iy = imfilter(double(I), hy, 'replicate');
    Ix = imfilter(double(I), hx, 'replicate');
    gradmag = sqrt(Ix.^2 + Iy.^2);
    % figure, imshow(gradmag,[]), title('Gradient magnitude (gradmag)')
    gradmag2 = imimposemin(gradmag, bgm | fgm4);
    % final watershed
    L = watershed(gradmag2);

    % Result visulization #1
    I4 = I;
    I4(imdilate(L == 0, ones(3, 3)) | bgm | fgm4) = 255;
    % I4(imdilate(DL == 0, ones(3, 3)) | bgm) = 255;
   figure, imshow(I4)
   title('Markers and object boundaries superimposed on original image (I4)')

    % Result visulization #2
%     Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
%     figure, imshow(Lrgb), title('Colored watershed label matrix (Lrgb)')
%    imwrite(Lrgb, ['.\TrenchImages\Trench\',imgname,'_',int2str(PL),'b.jpg']);
    % Result visulization #3
    %figure, imshow(I), hold on
    %himage = imshow(Lrgb);
    %set(himage, 'AlphaData', 0.3);
    %title('Lrgb superimposed transparently on original image')
  %  imwrite(I, ['.\TrenchImages\Trench\',imgname,'_',int2str(PL),'a.jpg']);