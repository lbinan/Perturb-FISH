%this script generates small images of individal segmented cells with
%transcripts in red, GFP (indicating a perturbed cell) in cyan, and saves
%images in folders according to local density
recimage=zeros(size(myimage(:,:,1)));
for k=1:size(stats,1)
    for j=1:size(stats(k).PixelIdxList,1)
        recimage(stats(k).PixelIdxList(j))=1;
    end
end

 recimage=(imdilate(recimage,strel('disk',3))-recimage)>0;
for i=1:size(numberofneighbors,1)
    neighbors=numberofneighbors(i)
%     GFP=gfp(i)
%     TNF=tnf(i)
%     IL1a=il1a(i)
    whereimage=zeros(size(myimage(:,:,1)));
    
    for j=1:size(stats(i).PixelIdxList,1)
    whereimage(stats(i).PixelIdxList(j))=1;
    end
    whereimage=imbinarize(whereimage);
    minY=min(find(sum(whereimage,1)>0))-40;
    maxY=max(find(sum(whereimage,1)>0))+40;
    minX=min(find(sum(whereimage,2)>0))-40;
    maxX=max(find(sum(whereimage,2)>0))+40;
   %TNF
    r=myimage(max(1,minX):min(2960,maxX),max(1,minY):min(2960,maxY),4)/7000;%3 15000
    g=myimage(max(1,minX):min(2960,maxX),max(1,minY):min(2960,maxY),2)/2000;
    croprecimage=recimage(max(1,minX):min(2960,maxX),max(1,minY):min(2960,maxY));
    g=g-r;
    redimage=0*g+r+0.6510*croprecimage;
    blueimage=0.4392*g+0.8078*croprecimage;
    greenimage=0.520*g+0.8902*croprecimage;
    thisImage=cat(3,redimage, greenimage, blueimage);
    if neighbors<neighborlowthresh
       if lowtnf<502
        lowtnf=lowtnf+1;
    	imwrite(thisImage,strcat(path,'/tnf/lowdensity/im',num2str(i),'.jpeg'))
       end
    elseif hightnf<300 & neighbors>neighborhighthresh
        hightnf=hightnf+1;
        imwrite(thisImage,strcat(path,'/tnf/highdensity/im',num2str(i),'.jpeg'))
    end
    % il1a
    r=myimage(max(1,minX):min(2960,maxX),max(1,minY):min(2960,maxY),3)/12000;
    g=myimage(max(1,minX):min(2960,maxX),max(1,minY):min(2960,maxY),2)/2000;
    croprecimage=recimage(max(1,minX):min(2960,maxX),max(1,minY):min(2960,maxY));
    g=g-r;
    redimage=0*g+r+0.6510*croprecimage;
    blueimage=0.4392*g+0.8078*croprecimage;
    greenimage=0.520*g+0.8902*croprecimage;
    thisImage=cat(3,redimage, greenimage, blueimage);
    if neighbors<neighborlowthresh
        if lowil1a<500
            lowil1a=lowil1a+1;
            imwrite(thisImage,strcat(path,'/il1a/lowdensity/im',num2str(i),'.jpeg'))
        end
    elseif highil1a<300 & neighbors>neighborhighthresh
        highil1a=highil1a+1;
        imwrite(thisImage,strcat(path,'/il1a/highdensity/im',num2str(i),'.jpeg'))
    end
end