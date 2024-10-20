function spots=decodeThisBinaryImageTUMORHW2(bw)
flatBW=sum(bw,3)>1;figure, imshow(flatBW)
allspotsCentroid=regionprops(flatBW,'Centroid');%get coordinates of putative guide spots
keepspotsIdx=zeros(1,size(allspotsCentroid,1));
spots=[];
checkrepeats=[1,2,3];%deprecated, used to check for spots that
% consitently show up the the same channel, was fixed in the lab by
% changing the codebook
checkrepeats=[checkrepeats,checkrepeats,checkrepeats,checkrepeats,checkrepeats];
for thisspot=1:size(allspotsCentroid,1)
    y=uint16(allspotsCentroid(thisspot).Centroid(1));
    x=uint16(allspotsCentroid(thisspot).Centroid(2));
    if x<26 | x>(2048-48) | y<26 | y>(2048-48)%removing the outside of the region
        signal(thisspot)=0;
    else
        signal=sum(max(max(bw(x-8:x+8,y-8:y+8,:))));%calculating the hamming weight of putative spot
    end
    if signal>2 & signal <6 %if the spot has at most one bit error, we look for its identity
        keepspotsIdx(thisspot)=1;
        spotvalues=[];
        for i=1:15
            spotvalues=[spotvalues;max(max(bw(max(1,x-8):min(2048,x+8),max(1,y-8):min(2048,y+8),i)))];
        end
        if max(size(unique(spotvalues.*checkrepeats)))>2
            spots=[spots;[double(x),double(y),spotvalues']];
        end
    end
end
