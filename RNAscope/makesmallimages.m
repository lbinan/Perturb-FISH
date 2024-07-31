%reads all raw images from RNA scope experiements for each gene we
%validated
mypath='/broad/clearylab/Users/Loic/rnascope/rnascope/nfkbi1high10x/20230824_012126_364'
savepath='/broad/clearylab/Users/Loic/rnascope/countspotsbetter'   
il1ahighperturb=[];
il1alowperturb=[];
tnfhighperturb=[];
tnflowperturb=[];
il1ahighcontrol=[];
il1alowcontrol=[];
tnfhighcontrol=[];
tnflowcontrol=[];
summary=[]
lowil1a=0;
highil1a=0;
lowtnf=0;
hightnf=0;
neighborlowthresh=2;
neighborhighthresh=3;

for jj=1:25
   	if jj~=16 & jj~=15 & jj~=21 %don't analyse out of focus images
    myimage=[];
    for i=1:4
        if jj<11
            myimage(:,:,i)=imread(fullfile(mypath,strcat('Point000',num2str(jj-1),'_ChannelSBS_DAPI,FITC-Penta,Alexa647,SS_Cy3_Seq000',num2str(jj-1),'.tiff')),i);
        else
            myimage(:,:,i)=imread(fullfile(mypath,strcat('Point00',num2str(jj-1),'_ChannelSBS_DAPI,FITC-Penta,Alexa647,SS_Cy3_Seq00',num2str(jj-1),'.tiff')),i);
        end
    end
    bw=imbinarize(myimage(:,:,1),7000);
    stats=regionprops(bw,'Centroid','PixelIdxList');
    seeds=zeros(size(bw));
    for i=1:size(stats,1)
        seeds(floor(stats(i).Centroid(2)),floor(stats(i).Centroid(1)))=1;
    end
    seeds=imdilate(seeds,strel('disk',5));
    bw=imbinarize(myimage(:,:,1)/5000,'adaptive');
    bw=bwareaopen(bw,100);
    D = bwdist(~bw);
    D = -D;
    gmag2 = imimposemin(D, seeds);
    mymask=watershed(gmag2);
    mymask(~bw)=0;
    toobig=bwareaopen(mymask>0,2000)|bwpropfilt(mymask>0,'Solidity',[0 0.88]);%we identify objects too big to be cells and split them again
    mymask(toobig)=0;
    D = bwdist(~toobig);
    D = -D;
    I3 = imhmin(D,0.5);
    resplitnuclei = watershed(I3);
    resplitnuclei(~toobig)=0;
    mymask=(mymask>0)+(resplitnuclei>0);
    mymask=bwareaopen(mymask,300);
    mymask=mymask-bwareaopen(mymask,3000);
    %% segment the whole cell
    stats=regionprops(mymask>0,'Centroid');
    seeds=zeros(size(bw));
    for i=1:size(stats,1)
        seeds(floor(stats(i).Centroid(2)),floor(stats(i).Centroid(1)))=1;
    end
    seeds=imdilate(seeds,strel('disk',10));
    cellarea=(myimage(:,:,3)+myimage(:,:,4))>2200;
    cellarea=bwareaopen(cellarea,500);%nake it 100?
    imshowpair(cellarea,seeds)
    D = bwdist(~cellarea);
    D = -D;
    gmag2 = imimposemin(D, seeds);

    mymask=watershed(gmag2);
    mymask(~cellarea)=0;
    %figure, imshowpair(mymask>0,seeds)
    gfp=myimage(:,:,2)>2000;
    %figure, imshow(gfp)
    gfp=bwareaopen(gfp,500);
    %figure, imshow(gfp)
    D = bwdist(~gfp);
    D = -D;
    gmag2 = imimposemin(D, seeds);

    perturbedmask=watershed(gmag2);
    perturbedmask(~gfp)=0;
    %figure, imshowpair(perturbedmask>0,mymask>0)
    perturbedmask=perturbedmask>0;
    mymask=mymask>0;
    mymask(imdilate(perturbedmask,strel('disk',2)))=0;
    mymask=mymask|perturbedmask;
    %figure, imshow(mymask)
    mymask=bwareaopen(mymask,500);
    %mymask=bwpropfilt(mymask,'EulerNumber',[0.8 1]); 
    mymask=bwpropfilt(mymask,'solidity',[0.7 1]);mymask=(mymask-bwareaopen(mymask,4000))>0;
    if jj==4 | jj==5 |jj==14 |jj==25
        mymask(:,1:800)=0;
    end
    if jj==22 | jj==23 | jj==24
        mymask(1:800,:)=0;
    end
    %%
    % %figure, imshowpair(toobig,mymask>0)
    % saturated=myimage(:,:,3)>50000|myimage(:,:,4)>50000;
    % saturated=bwareaopen(saturated,50);
    % saturated=imdilate(saturated,strel('disk',10));
    % satpixels=find(saturated);
    labeledimage=bwlabel(mymask);
    stats=regionprops(mymask>0,'Centroid','PixelIdxList');
    % gfp=zeros(size(stats,1),1);
    % il1a=zeros(size(stats,1),1);
    % tnf=zeros(size(stats,1),1);
    % perturbedcells=myimage(:,:,2);
    % il1aImage=myimage(:,:,3)/1200.*myimage(:,:,3)/1200;
    % il1aImage=il1aImage>26;%22
    % coor=pkfnd(imgaussfilt(myimage(:,:,3),1)-imgaussfilt(myimage(:,:,3),40),5000,2); 
    % il1aImage=zeros(2960,2960);
    % for j=1:size(coor,1)
    % il1aImage(coor(j,2),coor(j,1))=1;
    % end
    % tnfImage=myimage(:,:,4)/1200.*myimage(:,:,4)/1200;
    % tnfImage=tnfImage>8;%6

    % coor=pkfnd(imgaussfilt(myimage(:,:,4),1)-imgaussfilt(myimage(:,:,4),50),1500,2);
    % tnfImage=zeros(2960,2960);
    % for j=1:size(coor,1)
    % tnfImage(coor(j,2),coor(j,1))=1;
    % end
    % % saturatedcells=zeros(size(stats,1),1);
    % for i=1:size(gfp,1)
    % %     if size(find(ismember(satpixels,stats(i).PixelIdxList)),1)==0
    %     gfp(i)=mean(perturbedcells(stats(i).PixelIdxList));
    %     il1a(i)=sum(il1aImage(stats(i).PixelIdxList));
    %     tnf(i)=sum(tnfImage(stats(i).PixelIdxList));
    % %     else
    % %         saturatedcells(i)=1;
    % %     end
    % end
    % gfp=gfp(~logical(saturatedcells));
    % il1a=il1a(~logical(saturatedcells));
    % tnf=tnf(~logical(saturatedcells));

    %%figure,histogram(gfp)
    % maskperturbed=zeros(size(bw));
    % maskUNperturbed=zeros(size(bw));
    % for i=1:size(gfp,1)
    %     if size(find(ismember(satpixels,stats(i).PixelIdxList)),1)==0
    %     if gfp(i)>2000
    %         maskperturbed(stats(i).PixelIdxList)=1;
    %     elseif gfp(i)<500
    %         maskUNperturbed(stats(i).PixelIdxList)=1;
    %     end
    %     end
    % end
    % maskperturbed=maskperturbed>0;
    % maskUNperturbed=maskUNperturbed>0;
    %%figure, imshow(maskperturbed)
    %%figure, imshow(myimage(:,:,2),[])
    %%figure, imshowpair(mymask>0,maskperturbed)
    % perturbedindeces=gfp>1000;
    % unperturbedindeces=gfp<500;
    % 
    % perturbed=zeros(sum(perturbedindeces),2);
    % perturbed(:,1)=il1a(perturbedindeces);
    % perturbed(:,2)=tnf(perturbedindeces);
    % unperturbed=zeros(sum(unperturbedindeces),2);
    % unperturbed(:,1)=il1a(unperturbedindeces);
    % unperturbed(:,2)=tnf(unperturbedindeces);
    % mean(perturbed)
    % mean(unperturbed)

    % centroids=zeros(size(stats,1),2);
    % for i=1:size(stats,1)
    %     if size(find(ismember(satpixels,stats(i).PixelIdxList)),1)==0
    %     centroids(i,:)=stats(i).Centroid;
    %     end
    % end
    numberofneighbors=zeros(size(stats,1),1);
    parfor k=1:size(stats,1)%find neighbors by finding cells that overlap with morphologically dilated cell
    %   if size(find(ismember(satpixels,stats(k).PixelIdxList)),1)==0
        thiscellimage=zeros(size(myimage(:,:,1)));
        for j=1:size(stats(k).PixelIdxList,1)
            thiscellimage(stats(k).PixelIdxList(j))=1;
        end
        thiscellimage=imdilate(thiscellimage,strel('disk',10));
        numberofneighbors(k)=size(unique(labeledimage(thiscellimage>0)),1)-2;
    %     end
    end
    % numberofneighbors=numberofneighbors(~logical(saturatedcells));
    % distances=pdist2(centroids,centroids);
    % distances=distances<50&distances~=0;
    % distances(perturbedindeces,:)=0;
    % numberofneighbors=sum(distances);
    % cellsIDs=[1:size(distances,2)];
    % distances=distances.*cellsIDs;
    % numberofneighbors=numberofneighbors';
    % highperturb=(numberofneighbors>3)&perturbedindeces;
    % lowperturb=(numberofneighbors<2)&perturbedindeces;
    % highcontrol=(numberofneighbors>3)&unperturbedindeces; 
    % lowcontrol=(numberofneighbors<2)&unperturbedindeces;
    % il1ahighperturb=[il1ahighperturb;il1a(highperturb)];
    % il1alowperturb=[il1alowperturb;il1a(lowperturb)];
    % tnfhighperturb=[tnfhighperturb;tnf(highperturb)];
    % tnflowperturb=[tnflowperturb;tnf(lowperturb)];
    % 
    % il1ahighcontrol=[il1ahighcontrol;il1a(highcontrol)];
    % il1alowcontrol=[il1alowcontrol;il1a(lowcontrol)];
    % tnfhighcontrol=[tnfhighcontrol;tnf(highcontrol)];
    % tnflowcontrol=[tnflowcontrol;tnf(lowcontrol)];
    % summary=[summary;[il1a,tnf,gfp,numberofneighbors,jj*ones(size(numberofneighbors,1),1)]];
    printimages%calls script that generates miniature images of each cell
    end
end

% writematrix(summary,fullfile(savepath,'countsfornfkb1_1.csv'))
%%
%repeat everything for each guide RNA (here, same gene, guide number 2, then
%sequenitally all genes
mypath='/broad/clearylab/Users/Loic/rnascope/rnascope/nfkbi2lhigh10x/20230824_012615_237'
il1ahighperturb=[];
il1alowperturb=[];
tnfhighperturb=[];
tnflowperturb=[];
il1ahighcontrol=[];
il1alowcontrol=[];
tnfhighcontrol=[];
tnflowcontrol=[];
summary=[];

for jj=1:25
     if jj~=4 & jj~=20 & jj~=1& jj~=3 & jj~=5 & jj~=7  & jj~=24 & jj~=25
myimage=[];
for i=1:4
    if jj<11
        myimage(:,:,i)=imread(fullfile(mypath,strcat('Point000',num2str(jj-1),'_ChannelSBS_DAPI,FITC-Penta,Alexa647,SS_Cy3_Seq000',num2str(jj-1),'.tiff')),i);
    else
        myimage(:,:,i)=imread(fullfile(mypath,strcat('Point00',num2str(jj-1),'_ChannelSBS_DAPI,FITC-Penta,Alexa647,SS_Cy3_Seq00',num2str(jj-1),'.tiff')),i);
    end
end
bw=imbinarize(myimage(:,:,1),7000);
%%figure, imshow(bw)

stats=regionprops(bw,'Centroid','PixelIdxList');
seeds=zeros(size(bw));
for i=1:size(stats,1)
    seeds(floor(stats(i).Centroid(2)),floor(stats(i).Centroid(1)))=1;
end
seeds=imdilate(seeds,strel('disk',5));
%%figure, imshow(seeds);

bw=imbinarize(myimage(:,:,1)/5000,'adaptive');
%%figure, imshow(bw)
bw=bwareaopen(bw,100);
%%figure, imshow(bw)
D = bwdist(~bw);
D = -D;
gmag2 = imimposemin(D, seeds);

mymask=watershed(gmag2);
mymask(~bw)=0;
%%figure, imshow(mymask>0)
%%figure, imshowpair(mymask>0,seeds)

toobig=bwareaopen(mymask>0,2000)|bwpropfilt(mymask>0,'Solidity',[0 0.88]);
mymask(toobig)=0;

D = bwdist(~toobig);
D = -D;
I3 = imhmin(D,0.5);
resplitnuclei = watershed(I3);
resplitnuclei(~toobig)=0;
% %figure, imshowpair(resplitnuclei>0,mymask>0)
mymask=(mymask>0)+(resplitnuclei>0);
mymask=bwareaopen(mymask,300);
mymask=mymask-bwareaopen(mymask,3000);
%% segment the whole cell
stats=regionprops(mymask>0,'Centroid');
seeds=zeros(size(bw));
for i=1:size(stats,1)
    seeds(floor(stats(i).Centroid(2)),floor(stats(i).Centroid(1)))=1;
end
seeds=imdilate(seeds,strel('disk',10));
cellarea=(myimage(:,:,3)+myimage(:,:,4))>2200;
cellarea=bwareaopen(cellarea,500);
imshowpair(cellarea,seeds)
D = bwdist(~cellarea);
D = -D;
gmag2 = imimposemin(D, seeds);

mymask=watershed(gmag2);
mymask(~cellarea)=0;
%figure, imshowpair(mymask>0,seeds)
gfp=myimage(:,:,2)>2000;
%figure, imshow(gfp)
gfp=bwareaopen(gfp,500);
%figure, imshow(gfp)
D = bwdist(~gfp);
D = -D;
gmag2 = imimposemin(D, seeds);

perturbedmask=watershed(gmag2);
perturbedmask(~gfp)=0;
%figure, imshowpair(perturbedmask>0,mymask>0)
perturbedmask=perturbedmask>0;
mymask=mymask>0;
mymask(imdilate(perturbedmask,strel('disk',2)))=0;
mymask=mymask|perturbedmask;
%figure, imshow(mymask)
mymask=bwareaopen(mymask,500);
%mymask=bwpropfilt(mymask,'EulerNumber',[0.8 1]); 
mymask=bwpropfilt(mymask,'solidity',[0.7 1]);mymask=(mymask-bwareaopen(mymask,4000))>0;
%%
% %figure, imshowpair(toobig,mymask>0)
% labeledimage=bwlabel(mymask);
% saturated=myimage(:,:,3)>50000|myimage(:,:,4)>50000;
% saturated=bwareaopen(saturated,50);
% saturated=imdilate(saturated,strel('disk',10));
% satpixels=find(saturated);
labeledimage=bwlabel(mymask);
stats=regionprops(mymask>0,'Centroid','PixelIdxList');
% gfp=zeros(size(stats,1),1);
% il1a=zeros(size(stats,1),1);
% tnf=zeros(size(stats,1),1);
% perturbedcells=myimage(:,:,2);
% il1aImage=myimage(:,:,3)/1200.*myimage(:,:,3)/1200;
% il1aImage=il1aImage>26;%22
% coor=pkfnd(imgaussfilt(myimage(:,:,3),1)-imgaussfilt(myimage(:,:,3),40),5000,2); 
% il1aImage=zeros(2960,2960);
% for j=1:size(coor,1)
% il1aImage(coor(j,2),coor(j,1))=1;
% end
% tnfImage=myimage(:,:,4)/1200.*myimage(:,:,4)/1200;
% tnfImage=tnfImage>8;%6

% coor=pkfnd(imgaussfilt(myimage(:,:,4),1)-imgaussfilt(myimage(:,:,4),50),1500,2);
% tnfImage=zeros(2960,2960);
% for j=1:size(coor,1)
% tnfImage(coor(j,2),coor(j,1))=1;
% end
% saturatedcells=zeros(size(stats,1),1);
% for i=1:size(gfp,1)
% %     if size(find(ismember(satpixels,stats(i).PixelIdxList)),1)==0
%     gfp(i)=mean(perturbedcells(stats(i).PixelIdxList));
%     il1a(i)=sum(il1aImage(stats(i).PixelIdxList));
%     tnf(i)=sum(tnfImage(stats(i).PixelIdxList));
% %     else
% %         saturatedcells(i)=1;
% %     end
% end
% gfp=gfp(~logical(saturatedcells));
% il1a=il1a(~logical(saturatedcells));
% tnf=tnf(~logical(saturatedcells));

%%figure,histogram(gfp)
% maskperturbed=zeros(size(bw));
% maskUNperturbed=zeros(size(bw));
% for i=1:size(gfp,1)
%     if size(find(ismember(satpixels,stats(i).PixelIdxList)),1)==0
%     if gfp(i)>2000
%         maskperturbed(stats(i).PixelIdxList)=1;
%     elseif gfp(i)<500
%         maskUNperturbed(stats(i).PixelIdxList)=1;
%     end
%     end
% end
% maskperturbed=maskperturbed>0;
% maskUNperturbed=maskUNperturbed>0;
%%figure, imshow(maskperturbed)
%%figure, imshow(myimage(:,:,2),[])
%%figure, imshowpair(mymask>0,maskperturbed)
% perturbedindeces=gfp>1000;
% unperturbedindeces=gfp<500;
% 
% perturbed=zeros(sum(perturbedindeces),2);
% perturbed(:,1)=il1a(perturbedindeces);
% perturbed(:,2)=tnf(perturbedindeces);
% unperturbed=zeros(sum(unperturbedindeces),2);
% unperturbed(:,1)=il1a(unperturbedindeces);
% unperturbed(:,2)=tnf(unperturbedindeces);
% mean(perturbed)
% mean(unperturbed)

% centroids=zeros(size(stats,1),2);
% for i=1:size(stats,1)
%     if size(find(ismember(satpixels,stats(i).PixelIdxList)),1)==0
%     centroids(i,:)=stats(i).Centroid;
%     end
% end
numberofneighbors=zeros(size(stats,1),1);
parfor k=1:size(stats,1)
%     if size(find(ismember(satpixels,stats(k).PixelIdxList)),1)==0
    thiscellimage=zeros(size(myimage(:,:,1)));
    for j=1:size(stats(k).PixelIdxList,1)
    thiscellimage(stats(k).PixelIdxList(j))=1;
    end
    thiscellimage=imdilate(thiscellimage,strel('disk',10));
    numberofneighbors(k)=size(unique(labeledimage(thiscellimage>0)),1)-2;
%     end
end
% numberofneighbors=numberofneighbors(~logical(saturatedcells));
% distances=pdist2(centroids,centroids);
% distances=distances<50&distances~=0;
% distances(perturbedindeces,:)=0;
% numberofneighbors=sum(distances);
% cellsIDs=[1:size(distances,2)];
% distances=distances.*cellsIDs;
% numberofneighbors=numberofneighbors';
% highperturb=(numberofneighbors>3)&perturbedindeces;
% lowperturb=(numberofneighbors<2)&perturbedindeces;
% highcontrol=(numberofneighbors>3)&unperturbedindeces; 
% lowcontrol=(numberofneighbors<2)&unperturbedindeces;
% il1ahighperturb=[il1ahighperturb;il1a(highperturb)];
% il1alowperturb=[il1alowperturb;il1a(lowperturb)];
% tnfhighperturb=[tnfhighperturb;tnf(highperturb)];
% tnflowperturb=[tnflowperturb;tnf(lowperturb)];
% 
% il1ahighcontrol=[il1ahighcontrol;il1a(highcontrol)];
% il1alowcontrol=[il1alowcontrol;il1a(lowcontrol)];
% tnfhighcontrol=[tnfhighcontrol;tnf(highcontrol)];
% tnflowcontrol=[tnflowcontrol;tnf(lowcontrol)];
% summary=[summary;[il1a,tnf,gfp,numberofneighbors,jj*ones(size(numberofneighbors,1),1)]]; 
printimages
     end
end
% writematrix(summary,fullfile(savepath,'countsfornfkb1_2.csv'))
