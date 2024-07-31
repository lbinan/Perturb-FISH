calciumclusterallidxitcalciumZscoreONEatAtime;
idxcalciumwell1=calciumclusterallidxitcalciumZscoreONEatAtime(calciumclusterallidxitcalciumZscoreONEatAtime(:,1)==21,:);
pooledmerfish=[];
idx=ismember(merfishcounttableS2region0rightIndex(:,1),idxcalciumwell1(:,2));
pooledmerfish=[pooledmerfish;merfishcounttableS2region0rightIndex(idx,:)];

idxcalciumwell2=calciumclusterallidxitcalciumZscoreONEatAtime(calciumclusterallidxitcalciumZscoreONEatAtime(:,1)==22,:);
idx=ismember(merfishcounttableS2region1rightIndex(:,1),idxcalciumwell2(:,2));
pooledmerfish=[pooledmerfish;merfishcounttableS2region1rightIndex(idx,:)];

idxcalciumwell3=calciumclusterallidxitcalciumZscoreONEatAtime(calciumclusterallidxitcalciumZscoreONEatAtime(:,1)==31,:);
idx=ismember(merfishcounttableS3region0rightIndex(:,1),idxcalciumwell3(:,2));
pooledmerfish=[pooledmerfish;merfishcounttableS3region0rightIndex(idx,:)];

idxcalciumwell4=calciumclusterallidxitcalciumZscoreONEatAtime(calciumclusterallidxitcalciumZscoreONEatAtime(:,1)==32,:);
idx=ismember(merfishcounttableS3region1rightindex(:,1),idxcalciumwell4(:,2));
pooledmerfish=[pooledmerfish;merfishcounttableS3region1rightindex(idx,:)];

idxcalciumwell5=calciumclusterallidxitcalciumZscoreONEatAtime(calciumclusterallidxitcalciumZscoreONEatAtime(:,1)==4,:);
idx=ismember(merfishcounttableS4rightIndex(:,1),idxcalciumwell5(:,2));
pooledmerfish=[pooledmerfish;merfishcounttableS4rightIndex(idx,:)];
writematrix(pooledmerfish,'\\helium\broad_clearylab\Users\Loic\AsdRevstables\pooledmerfishmatchedtocalcium.csv')

pooledmerfish=pooledmerfish(:,3:end);
totalcounts=sum(pooledmerfish,2);
badcells=totalcounts<200;
badcells=badcells|(totalcounts>1900);
pooledmerfish=pooledmerfish(~badcells,:);
calcium=calciumclusterallidxitcalciumZscoreONEatAtime(~badcells,:);
pooledmerfish=pooledmerfish(:,1:485);
zscored=pooledmerfish;
zscored(calcium(:,1)==21,:)=(pooledmerfish(calcium(:,1)==21,:)-mean(pooledmerfish(calcium(:,1)==21,:),1))./std(pooledmerfish(calcium(:,1)==21,:),1);
zscored(calcium(:,1)==22,:)=(pooledmerfish(calcium(:,1)==22,:)-mean(pooledmerfish(calcium(:,1)==22,:),1))./std(pooledmerfish(calcium(:,1)==22,:),1);
zscored(calcium(:,1)==31,:)=(pooledmerfish(calcium(:,1)==31,:)-mean(pooledmerfish(calcium(:,1)==31,:),1))./std(pooledmerfish(calcium(:,1)==31,:),1);
zscored(calcium(:,1)==32,:)=(pooledmerfish(calcium(:,1)==32,:)-mean(pooledmerfish(calcium(:,1)==32,:),1))./std(pooledmerfish(calcium(:,1)==32,:),1);
zscored(calcium(:,1)==4,:)=(pooledmerfish(calcium(:,1)==4,:)-mean(pooledmerfish(calcium(:,1)==4,:),1))./std(pooledmerfish(calcium(:,1)==4,:),1);

for j=1:max(calcium(:,3))
    Zscoredgenesignature(j,:)=mean(zscored(calcium(:,3)==j,:));
    controlsignature(j,:)=median(zscored(calcium(:,3)~=j,:));
end
% figure, heatmap(Zscoredgenesignature,'Colormap',redblue(256))
expressedgenesID=mean(pooledmerfish,1)>1;
figure, heatmap(Zscoredgenesignature(:,expressedgenesID),'Colormap',redblue(256))
clusternames=["large peak","inactive","large early transient","small peak","step","oscillating","step2","delayed"];
figure, heatmap(geneNames(expressedgenesID),clusternames,Zscoredgenesignature(:,expressedgenesID),'Colormap',redblue(256))

merfishedcalcium=calcium;
merfishedcalcium(merfishedcalcium(:,3)==7,3)=5;
merfishedcalcium(merfishedcalcium(:,3)==6,3)=2;
merfishedcalcium(merfishedcalcium(:,3)==8,3)=6;
Zscoredgenesignature=[];
for j=1:max(merfishedcalcium(:,3))
    Zscoredgenesignature(j,:)=mean(zscored(merfishedcalcium(:,3)==j,:));
    
end
% figure, heatmap(Zscoredgenesignature,'Colormap',redblue(256))
clusternames=["large peak","inactive","large early transient","small peak","step","delayed"];
figure, heatmap(geneNames(expressedgenesID),clusternames,Zscoredgenesignature(:,expressedgenesID),'Colormap',redblue(256))

for i=1:485
    p(i)=anova1(zscored(:,i),calcium(:,3),'off');
end
sig=p<0.05;
figure, heatmap(geneNames(expressedgenesID&sig),clusternames,Zscoredgenesignature(:,expressedgenesID&sig),'Colormap',redblue(256))
figure, heatmap(geneNames(sig),clusternames,Zscoredgenesignature(:,sig),'Colormap',redblue(256))
figure, heatmap(geneNames(sig&max(abs(Zscoredgenesignature))>0.1),clusternames,Zscoredgenesignature(:,sig&max(abs(Zscoredgenesignature))>0.1),'Colormap',redblue(256))
figure, heatmap(geneNames(strcmp(geneNames,"AXL")),clusternames,Zscoredgenesignature(:,strcmp(geneNames,"AXL")),'Colormap',redblue(256))
figure, heatmap(geneNames(strcmp(geneNames,"HDLBP")),clusternames,Zscoredgenesignature(:,strcmp(geneNames,"HDLBP")),'Colormap',redblue(256))

%% norm to total counts
merfishedcalcium=calcium;
merfishedcalcium(merfishedcalcium(:,3)==7,3)=5;
merfishedcalcium(merfishedcalcium(:,3)==6,3)=2;
merfishedcalcium(merfishedcalcium(:,3)==8,3)=6;
zscored=pooledmerfish./sum(pooledmerfish,2);
zscored(calcium(:,1)==21,:)=(zscored(calcium(:,1)==21,:)-mean(zscored(calcium(:,1)==21,:),1))./std(zscored(calcium(:,1)==21,:),1);
zscored(calcium(:,1)==22,:)=(zscored(calcium(:,1)==22,:)-mean(zscored(calcium(:,1)==22,:),1))./std(zscored(calcium(:,1)==22,:),1);
zscored(calcium(:,1)==31,:)=(zscored(calcium(:,1)==31,:)-mean(zscored(calcium(:,1)==31,:),1))./std(zscored(calcium(:,1)==31,:),1);
zscored(calcium(:,1)==32,:)=(zscored(calcium(:,1)==32,:)-mean(zscored(calcium(:,1)==32,:),1))./std(zscored(calcium(:,1)==32,:),1);
zscored(calcium(:,1)==4,:)=(zscored(calcium(:,1)==4,:)-mean(zscored(calcium(:,1)==4,:),1))./std(zscored(calcium(:,1)==4,:),1);


Zscoredgenesignature=[];
for j=1:max(merfishedcalcium(:,3))
    Zscoredgenesignature(j,:)=mean(zscored(merfishedcalcium(:,3)==j,:));
end
% figure, heatmap(Zscoredgenesignature,'Colormap',redblue(256))
clusternames=["large peak","inactive","large early transient","small peak","step","delayed"];
figure, heatmap(geneNames(expressedgenesID),clusternames,Zscoredgenesignature(:,expressedgenesID),'Colormap',redblue(256))

for i=1:485
    p(i)=anova1(zscored(:,i),calcium(:,3),'off');
end
sig=p<0.05;
figure, heatmap(geneNames(expressedgenesID&sig),clusternames,Zscoredgenesignature(:,expressedgenesID&sig),'Colormap',redblue(256))
figure, heatmap(geneNames(sig),clusternames,Zscoredgenesignature(:,sig),'Colormap',redblue(256))
title("p<0.05 , normlized to total counts per cell")
figure, heatmap(geneNames(sig&max(abs(Zscoredgenesignature))>0.1),clusternames,Zscoredgenesignature(:,sig&max(abs(Zscoredgenesignature))>0.1),'Colormap',redblue(256))
title("p<0.05 and |zscore|>0.1, normlized to total counts per cell")
figure, heatmap(geneNames(strcmp(geneNames,"AXL")),clusternames,Zscoredgenesignature(:,strcmp(geneNames,"AXL")),'Colormap',redblue(256))
figure, heatmap(geneNames(strcmp(geneNames,"HDLBP")),clusternames,Zscoredgenesignature(:,strcmp(geneNames,"HDLBP")),'Colormap',redblue(256))

totalcounts=sum(pooledmerfish,2);
for i=1:max(merfishedcalcium(:,3))
    totalexp(i)=mean(totalcounts(calcium(:,3)==i));
end
figure, heatmap("total counts",clusternames,totalexp','Colormap',redblue(256))


totalzscored=sum(zscored,2);
for i=1:max(merfishedcalcium(:,3))
    totalexpcorr(i)=mean(totalzscored(calcium(:,3)==i));
end
figure, heatmap("total counts",clusternames,totalexpcorr','Colormap',redblue(256))



ZscoredgenesignatureS1=[];
for j=1:max(merfishedcalcium(:,3))
    ZscoredgenesignatureS1(j,:)=mean(zscored(merfishedcalcium(:,3)==j&merfishedcalcium(:,1)==21,:));
end
ZscoredgenesignatureS2=[];
for j=1:max(merfishedcalcium(:,3))
    ZscoredgenesignatureS2(j,:)=mean(zscored(merfishedcalcium(:,3)==j&merfishedcalcium(:,1)==22,:));
end
ZscoredgenesignatureS3=[];
for j=1:max(merfishedcalcium(:,3))
    ZscoredgenesignatureS3(j,:)=mean(zscored(merfishedcalcium(:,3)==j&merfishedcalcium(:,1)==31,:));
end
ZscoredgenesignatureS4=[];
for j=1:max(merfishedcalcium(:,3))
    ZscoredgenesignatureS4(j,:)=mean(zscored(merfishedcalcium(:,3)==j&merfishedcalcium(:,1)==32,:));
end
ZscoredgenesignatureS5=[];
for j=1:max(merfishedcalcium(:,3))
    ZscoredgenesignatureS5(j,:)=mean(zscored(merfishedcalcium(:,3)==j&merfishedcalcium(:,1)==4,:));
end
signaturepersample=[reshape(ZscoredgenesignatureS1,[6*485,1]),reshape(ZscoredgenesignatureS2,[6*485,1]),reshape(ZscoredgenesignatureS3,[6*485,1]),reshape(ZscoredgenesignatureS4,[6*485,1]),reshape(ZscoredgenesignatureS5,[6*485,1])];
figure, 
subplot(2,2,1)
scatter(signaturepersample(:,1),signaturepersample(:,2),'.')
title("zscored expression per cluster in well 1 vs well 2")
subplot(2,2,2)
 scatter(signaturepersample(:,1),signaturepersample(:,3),'.')
title("zscored expression per cluster in well 1 vs well 3")
subplot(2,2,3)
 scatter(signaturepersample(:,1),signaturepersample(:,4),'.')
title("zscored expression per cluster in well 1 vs well 4")
subplot(2,2,4)
 scatter(signaturepersample(:,1),signaturepersample(:,5),'.')
title("zscored expression per cluster in well 1 vs well 5")


figure, 
subplot(2,2,1)
scatter(signaturepersample(:,2),signaturepersample(:,1),'.')
title("zscored expression per cluster in well 2 vs well 1")
subplot(2,2,2)
 scatter(signaturepersample(:,2),signaturepersample(:,3),'.')
title("zscored expression per cluster in well 2 vs well 3")
subplot(2,2,3)
 scatter(signaturepersample(:,2),signaturepersample(:,4),'.')
title("zscored expression per cluster in well 2 vs well 4")
subplot(2,2,4)
 scatter(signaturepersample(:,2),signaturepersample(:,5),'.')
title("zscored expression per cluster in well 2 vs well 5")

%% again not normalized to total counts per cell
zscored=pooledmerfish;
zscored(calcium(:,1)==21,:)=(zscored(calcium(:,1)==21,:)-mean(zscored(calcium(:,1)==21,:),1))./std(zscored(calcium(:,1)==21,:),1);
zscored(calcium(:,1)==22,:)=(zscored(calcium(:,1)==22,:)-mean(zscored(calcium(:,1)==22,:),1))./std(zscored(calcium(:,1)==22,:),1);
zscored(calcium(:,1)==31,:)=(zscored(calcium(:,1)==31,:)-mean(zscored(calcium(:,1)==31,:),1))./std(zscored(calcium(:,1)==31,:),1);
zscored(calcium(:,1)==32,:)=(zscored(calcium(:,1)==32,:)-mean(zscored(calcium(:,1)==32,:),1))./std(zscored(calcium(:,1)==32,:),1);
zscored(calcium(:,1)==4,:)=(zscored(calcium(:,1)==4,:)-mean(zscored(calcium(:,1)==4,:),1))./std(zscored(calcium(:,1)==4,:),1);


ZscoredgenesignatureS1=[];
for j=1:max(merfishedcalcium(:,3))
    ZscoredgenesignatureS1(j,:)=mean(zscored(merfishedcalcium(:,3)==j&merfishedcalcium(:,1)==21,:));
end
ZscoredgenesignatureS2=[];
for j=1:max(merfishedcalcium(:,3))
    ZscoredgenesignatureS2(j,:)=mean(zscored(merfishedcalcium(:,3)==j&merfishedcalcium(:,1)==22,:));
end
ZscoredgenesignatureS3=[];
for j=1:max(merfishedcalcium(:,3))
    ZscoredgenesignatureS3(j,:)=mean(zscored(merfishedcalcium(:,3)==j&merfishedcalcium(:,1)==31,:));
end
ZscoredgenesignatureS4=[];
for j=1:max(merfishedcalcium(:,3))
    ZscoredgenesignatureS4(j,:)=mean(zscored(merfishedcalcium(:,3)==j&merfishedcalcium(:,1)==32,:));
end
ZscoredgenesignatureS5=[];
for j=1:max(merfishedcalcium(:,3))
    ZscoredgenesignatureS5(j,:)=mean(zscored(merfishedcalcium(:,3)==j&merfishedcalcium(:,1)==4,:));
end
signaturepersample=[reshape(ZscoredgenesignatureS1,[6*485,1]),reshape(ZscoredgenesignatureS2,[6*485,1]),reshape(ZscoredgenesignatureS3,[6*485,1]),reshape(ZscoredgenesignatureS4,[6*485,1]),reshape(ZscoredgenesignatureS5,[6*485,1])];
figure, 
subplot(2,2,1)
scatter(signaturepersample(:,1),signaturepersample(:,2),'.')
title("zscored expression per cluster in well 1 vs well 2")
subplot(2,2,2)
 scatter(signaturepersample(:,1),signaturepersample(:,3),'.')
title("zscored expression per cluster in well 1 vs well 3")
subplot(2,2,3)
 scatter(signaturepersample(:,1),signaturepersample(:,4),'.')
title("zscored expression per cluster in well 1 vs well 4")
subplot(2,2,4)
 scatter(signaturepersample(:,1),signaturepersample(:,5),'.')
title("zscored expression per cluster in well 1 vs well 5")


figure, 
subplot(2,2,1)
scatter(signaturepersample(:,2),signaturepersample(:,1),'.')
title("zscored expression per cluster in well 2 vs well 1")
subplot(2,2,2)
 scatter(signaturepersample(:,2),signaturepersample(:,3),'.')
title("zscored expression per cluster in well 2 vs well 3")
subplot(2,2,3)
 scatter(signaturepersample(:,2),signaturepersample(:,4),'.')
title("zscored expression per cluster in well 2 vs well 4")
subplot(2,2,4)
 scatter(signaturepersample(:,2),signaturepersample(:,5),'.')
title("zscored expression per cluster in well 2 vs well 5")


figure, 
subplot(2,2,1)
scatter(ZscoredgenesignatureS1(1,:),ZscoredgenesignatureS2(1,:),'.')
title("zscored expression for cluster large peak in well 1 vs well 2")
subplot(2,2,2)
 scatter(ZscoredgenesignatureS1(1,:),ZscoredgenesignatureS2(1,:),'.')
title("zscored expression for cluster large peak in well 1 vs well 2")
subplot(2,2,3)
 scatter(ZscoredgenesignatureS1(1,:),ZscoredgenesignatureS4(1,:),'.')
title("zscored expression for cluster large peak in well 1 vs well 4")
subplot(2,2,4)
 scatter(ZscoredgenesignatureS1(1,:),ZscoredgenesignatureS5(1,:),'.')
title("zscored expression for cluster large peak in well 1 vs well 5")

figure, 
subplot(2,2,1)
scatter(ZscoredgenesignatureS1(2,:),ZscoredgenesignatureS2(2,:),'.')
title("zscored expression for cluster inactive in well 1 vs well 2")
subplot(2,2,2)
 scatter(ZscoredgenesignatureS1(2,:),ZscoredgenesignatureS2(2,:),'.')
title("zscored expression for cluster inactive in well 1 vs well 2")
subplot(2,2,3)
 scatter(ZscoredgenesignatureS1(2,:),ZscoredgenesignatureS4(2,:),'.')
title("zscored expression for cluster inactive in well 1 vs well 4")
subplot(2,2,4)
 scatter(ZscoredgenesignatureS1(2,:),ZscoredgenesignatureS5(2,:),'.')
title("zscored expression for cluster inactive in well 1 vs well 5")

figure, 
subplot(2,2,1)
scatter(ZscoredgenesignatureS1(3,:),ZscoredgenesignatureS2(3,:),'.')
title("zscored expression for cluster large early transient in well 1 vs well 2")
subplot(2,2,2)
 scatter(ZscoredgenesignatureS1(3,:),ZscoredgenesignatureS3(3,:),'.')
title("zscored expression for cluster llarge early transient in well 1 vs well 3")
subplot(2,2,3)
 scatter(ZscoredgenesignatureS1(3,:),ZscoredgenesignatureS4(3,:),'.')
title("zscored expression for cluster large early transient in well 1 vs well 4")
subplot(2,2,4)
 scatter(ZscoredgenesignatureS1(3,:),ZscoredgenesignatureS5(3,:),'.')
title("zscored expression for cluster large early transient in well 1 vs well 5")

figure, 
subplot(2,2,1)
scatter(ZscoredgenesignatureS1(4,:),ZscoredgenesignatureS2(4,:),'.')
title("zscored expression for cluster small peak in well 1 vs well 2")
subplot(2,2,2)
 scatter(ZscoredgenesignatureS1(4,:),ZscoredgenesignatureS3(4,:),'.')
title("zscored expression for cluster small peak in well 1 vs well 3")
subplot(2,2,3)
 scatter(ZscoredgenesignatureS1(4,:),ZscoredgenesignatureS4(4,:),'.')
title("zscored expression for cluster small peak in well 1 vs well 4")
subplot(2,2,4)
 scatter(ZscoredgenesignatureS1(4,:),ZscoredgenesignatureS5(4,:),'.')
title("zscored expression for cluster small peak in well 1 vs well 5")


figure, 
subplot(2,2,1)
scatter(ZscoredgenesignatureS1(5,:),ZscoredgenesignatureS2(5,:),'.')
title("zscored expression for cluster step in well 1 vs well 2")
subplot(2,2,2)
 scatter(ZscoredgenesignatureS1(5,:),ZscoredgenesignatureS3(5,:),'.')
title("zscored expression for cluster step in well 1 vs well 3")
subplot(2,2,3)
 scatter(ZscoredgenesignatureS1(5,:),ZscoredgenesignatureS4(5,:),'.')
title("zscored expression for cluster step in well 1 vs well 4")
subplot(2,2,4)
 scatter(ZscoredgenesignatureS1(5,:),ZscoredgenesignatureS5(5,:),'.')
title("zscored expression for cluster step in well 1 vs well 5")

figure, 
subplot(2,2,1)
scatter(ZscoredgenesignatureS1(6,:),ZscoredgenesignatureS2(6,:),'.')
title("zscored expression for cluster delayed in well 1 vs well 2")
subplot(2,2,2)
 scatter(ZscoredgenesignatureS1(6,:),ZscoredgenesignatureS3(6,:),'.')
title("zscored expression for cluster delayed in well 1 vs well 3")
subplot(2,2,3)
 scatter(ZscoredgenesignatureS1(6,:),ZscoredgenesignatureS4(6,:),'.')
title("zscored expression for cluster delayed in well 1 vs well 4")
subplot(2,2,4)
 scatter(ZscoredgenesignatureS1(6,:),ZscoredgenesignatureS5(6,:),'.')
title("zscored expression for cluster delayed in well 1 vs well 5")

%% repeat with normalized to toal counts
figure, 
subplot(2,2,1)
scatter(ZscoredgenesignatureS1(1,:),ZscoredgenesignatureS2(1,:),'.')
c=corrcoef(ZscoredgenesignatureS1(1,:),ZscoredgenesignatureS2(1,:));
d=corrcoef(ZscoredgenesignatureS1(1,abs(ZscoredgenesignatureS1(1,:))>0.07),ZscoredgenesignatureS2(1,abs(ZscoredgenesignatureS1(1,:))>0.07));
title(strcat("well 1 vs well 2, corr1 ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,2)
scatter(ZscoredgenesignatureS1(1,:),ZscoredgenesignatureS3(1,:),'.')
c=corrcoef(ZscoredgenesignatureS1(1,:),ZscoredgenesignatureS3(1,:));
d=corrcoef(ZscoredgenesignatureS1(1,abs(ZscoredgenesignatureS1(1,:))>0.07),ZscoredgenesignatureS3(1,abs(ZscoredgenesignatureS1(1,:))>0.07));
title(strcat("well 1 vs well 3, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,3)
 scatter(ZscoredgenesignatureS1(1,:),ZscoredgenesignatureS4(1,:),'.')
c=corrcoef(ZscoredgenesignatureS1(1,:),ZscoredgenesignatureS4(1,:));
d=corrcoef(ZscoredgenesignatureS1(1,abs(ZscoredgenesignatureS1(1,:))>0.07),ZscoredgenesignatureS4(1,abs(ZscoredgenesignatureS1(1,:))>0.07));
title(strcat("well 1 vs well 4, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,4)
 scatter(ZscoredgenesignatureS1(1,:),ZscoredgenesignatureS5(1,:),'.')
c=corrcoef(ZscoredgenesignatureS1(1,:),ZscoredgenesignatureS5(1,:));
d=corrcoef(ZscoredgenesignatureS1(1,abs(ZscoredgenesignatureS1(1,:))>0.07),ZscoredgenesignatureS5(1,abs(ZscoredgenesignatureS1(1,:))>0.07));
title(strcat("well 1 vs well 5, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
sgtitle('zscored expression for cluster large peak, normed to total counts') 
figure, 
subplot(2,2,1)
scatter(ZscoredgenesignatureS1(2,:),ZscoredgenesignatureS2(2,:),'.')
c=corrcoef(ZscoredgenesignatureS1(2,:),ZscoredgenesignatureS2(2,:));
d=corrcoef(ZscoredgenesignatureS1(2,abs(ZscoredgenesignatureS1(2,:))>0.07),ZscoredgenesignatureS2(2,abs(ZscoredgenesignatureS1(2,:))>0.07));
title(strcat("well 1 vs well 2, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,2)
scatter(ZscoredgenesignatureS1(2,:),ZscoredgenesignatureS3(2,:),'.')
c=corrcoef(ZscoredgenesignatureS1(2,:),ZscoredgenesignatureS3(2,:));
d=corrcoef(ZscoredgenesignatureS1(2,abs(ZscoredgenesignatureS1(2,:))>0.07),ZscoredgenesignatureS3(2,abs(ZscoredgenesignatureS1(2,:))>0.07));
title(strcat("well 1 vs well 3, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,3)
 scatter(ZscoredgenesignatureS1(2,:),ZscoredgenesignatureS4(2,:),'.')
c=corrcoef(ZscoredgenesignatureS1(2,:),ZscoredgenesignatureS4(2,:));
d=corrcoef(ZscoredgenesignatureS1(2,abs(ZscoredgenesignatureS1(2,:))>0.07),ZscoredgenesignatureS4(2,abs(ZscoredgenesignatureS1(2,:))>0.07));
title(strcat("well 1 vs well 4, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,4)
 scatter(ZscoredgenesignatureS1(2,:),ZscoredgenesignatureS5(2,:),'.')
c=corrcoef(ZscoredgenesignatureS1(2,:),ZscoredgenesignatureS5(2,:));
d=corrcoef(ZscoredgenesignatureS1(2,abs(ZscoredgenesignatureS1(2,:))>0.07),ZscoredgenesignatureS5(2,abs(ZscoredgenesignatureS1(2,:))>0.07));
title(strcat("well 1 vs well 5, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
sgtitle('zscored expression for cluster inactive cells, normed to total counts') 
figure, 
subplot(2,2,1)
scatter(ZscoredgenesignatureS1(3,:),ZscoredgenesignatureS2(3,:),'.')
c=corrcoef(ZscoredgenesignatureS1(3,:),ZscoredgenesignatureS2(3,:));
d=corrcoef(ZscoredgenesignatureS1(3,abs(ZscoredgenesignatureS1(3,:))>0.07),ZscoredgenesignatureS2(3,abs(ZscoredgenesignatureS1(3,:))>0.07));
title(strcat("well 1 vs well 2, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,2)
scatter(ZscoredgenesignatureS1(3,:),ZscoredgenesignatureS3(3,:),'.')
c=corrcoef(ZscoredgenesignatureS1(3,:),ZscoredgenesignatureS3(3,:));
d=corrcoef(ZscoredgenesignatureS1(3,abs(ZscoredgenesignatureS1(3,:))>0.07),ZscoredgenesignatureS3(3,abs(ZscoredgenesignatureS1(3,:))>0.07));
title(strcat("well 1 vs well 3, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,3)
 scatter(ZscoredgenesignatureS1(3,:),ZscoredgenesignatureS4(3,:),'.')
c=corrcoef(ZscoredgenesignatureS1(3,:),ZscoredgenesignatureS4(3,:));
d=corrcoef(ZscoredgenesignatureS1(3,abs(ZscoredgenesignatureS1(3,:))>0.07),ZscoredgenesignatureS4(3,abs(ZscoredgenesignatureS1(3,:))>0.07));
title(strcat("well 1 vs well 4, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,4)
 scatter(ZscoredgenesignatureS1(3,:),ZscoredgenesignatureS5(3,:),'.')
c=corrcoef(ZscoredgenesignatureS1(3,:),ZscoredgenesignatureS5(3,:));
d=corrcoef(ZscoredgenesignatureS1(3,abs(ZscoredgenesignatureS1(3,:))>0.07),ZscoredgenesignatureS5(3,abs(ZscoredgenesignatureS1(3,:))>0.07));
title(strcat("well 1 vs well 5, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
sgtitle('zscored expression for cluster large early transients, normed to total counts') 
figure, 
subplot(2,2,1)
scatter(ZscoredgenesignatureS1(4,:),ZscoredgenesignatureS2(4,:),'.')
c=corrcoef(ZscoredgenesignatureS1(4,:),ZscoredgenesignatureS2(4,:));
d=corrcoef(ZscoredgenesignatureS1(4,abs(ZscoredgenesignatureS1(4,:))>0.07),ZscoredgenesignatureS2(4,abs(ZscoredgenesignatureS1(4,:))>0.07));
title(strcat("well 1 vs well 2, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,2)
scatter(ZscoredgenesignatureS1(4,:),ZscoredgenesignatureS3(4,:),'.')
c=corrcoef(ZscoredgenesignatureS1(4,:),ZscoredgenesignatureS3(4,:));
d=corrcoef(ZscoredgenesignatureS1(4,abs(ZscoredgenesignatureS1(4,:))>0.07),ZscoredgenesignatureS3(4,abs(ZscoredgenesignatureS1(4,:))>0.07));
title(strcat("well 1 vs well 3, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,3)
 scatter(ZscoredgenesignatureS1(4,:),ZscoredgenesignatureS4(4,:),'.')
c=corrcoef(ZscoredgenesignatureS1(4,:),ZscoredgenesignatureS4(4,:));
d=corrcoef(ZscoredgenesignatureS1(4,abs(ZscoredgenesignatureS1(4,:))>0.07),ZscoredgenesignatureS4(4,abs(ZscoredgenesignatureS1(4,:))>0.07));
title(strcat("well 1 vs well 4, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,4)
 scatter(ZscoredgenesignatureS1(4,:),ZscoredgenesignatureS5(4,:),'.')
c=corrcoef(ZscoredgenesignatureS1(4,:),ZscoredgenesignatureS5(4,:));
d=corrcoef(ZscoredgenesignatureS1(4,abs(ZscoredgenesignatureS1(4,:))>0.07),ZscoredgenesignatureS2(4,abs(ZscoredgenesignatureS1(4,:))>0.07));
title(strcat("well 1 vs well 5, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
sgtitle('zscored expression for cluster small peak, normed to total counts') 


figure, 
subplot(2,2,1)
scatter(ZscoredgenesignatureS1(5,:),ZscoredgenesignatureS2(5,:),'.')
c=corrcoef(ZscoredgenesignatureS1(5,:),ZscoredgenesignatureS2(5,:));
d=corrcoef(ZscoredgenesignatureS1(5,abs(ZscoredgenesignatureS1(5,:))>0.07),ZscoredgenesignatureS2(5,abs(ZscoredgenesignatureS1(5,:))>0.07));
title(strcat("well 1 vs well 2, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,2)
scatter(ZscoredgenesignatureS1(5,:),ZscoredgenesignatureS3(5,:),'.')
c=corrcoef(ZscoredgenesignatureS1(5,:),ZscoredgenesignatureS3(5,:));
d=corrcoef(ZscoredgenesignatureS1(5,abs(ZscoredgenesignatureS1(5,:))>0.07),ZscoredgenesignatureS3(5,abs(ZscoredgenesignatureS1(5,:))>0.07));
title(strcat("well 1 vs well 3, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,3)
 scatter(ZscoredgenesignatureS1(5,:),ZscoredgenesignatureS4(5,:),'.')
c=corrcoef(ZscoredgenesignatureS1(5,:),ZscoredgenesignatureS4(5,:));
d=corrcoef(ZscoredgenesignatureS1(5,abs(ZscoredgenesignatureS1(5,:))>0.07),ZscoredgenesignatureS4(5,abs(ZscoredgenesignatureS1(5,:))>0.07));
title(strcat("well 1 vs well 4, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,4)
 scatter(ZscoredgenesignatureS1(5,:),ZscoredgenesignatureS5(5,:),'.')
c=corrcoef(ZscoredgenesignatureS1(5,:),ZscoredgenesignatureS5(5,:));
d=corrcoef(ZscoredgenesignatureS1(5,abs(ZscoredgenesignatureS1(5,:))>0.07),ZscoredgenesignatureS5(5,abs(ZscoredgenesignatureS1(5,:))>0.07));
title(strcat("well 1 vs well 5, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
sgtitle('zscored expression for cluster step, normed to total counts') 

figure, 
subplot(2,2,1)
scatter(ZscoredgenesignatureS1(6,:),ZscoredgenesignatureS2(6,:),'.')
c=corrcoef(ZscoredgenesignatureS1(6,:),ZscoredgenesignatureS2(6,:));
d=corrcoef(ZscoredgenesignatureS1(6,abs(ZscoredgenesignatureS1(6,:))>0.07),ZscoredgenesignatureS2(6,abs(ZscoredgenesignatureS1(6,:))>0.07));
title(strcat("well 1 vs well 2, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,2)
scatter(ZscoredgenesignatureS1(6,:),ZscoredgenesignatureS3(6,:),'.')
c=corrcoef(ZscoredgenesignatureS1(6,:),ZscoredgenesignatureS3(6,:));
d=corrcoef(ZscoredgenesignatureS1(6,abs(ZscoredgenesignatureS1(6,:))>0.07),ZscoredgenesignatureS3(6,abs(ZscoredgenesignatureS1(6,:))>0.07));
title(strcat("well 1 vs well 3, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,3)
 scatter(ZscoredgenesignatureS1(6,:),ZscoredgenesignatureS4(6,:),'.')
c=corrcoef(ZscoredgenesignatureS1(6,:),ZscoredgenesignatureS4(6,:));
d=corrcoef(ZscoredgenesignatureS1(6,abs(ZscoredgenesignatureS1(6,:))>0.07),ZscoredgenesignatureS4(6,abs(ZscoredgenesignatureS1(6,:))>0.07));
title(strcat("well 1 vs well 4, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,4)
 scatter(ZscoredgenesignatureS1(6,:),ZscoredgenesignatureS5(6,:),'.')
c=corrcoef(ZscoredgenesignatureS1(6,:),ZscoredgenesignatureS5(6,:));
d=corrcoef(ZscoredgenesignatureS1(6,abs(ZscoredgenesignatureS1(6,:))>0.07),ZscoredgenesignatureS5(6,abs(ZscoredgenesignatureS1(6,:))>0.07));
title(strcat("well 1 vs well 5, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
sgtitle('zscored expression for cluster delayed, normed to total counts') 

%% once more, no norm
figure, 
subplot(2,2,1)
scatter(ZscoredgenesignatureS1(1,:),ZscoredgenesignatureS2(1,:),'.')
c=corrcoef(ZscoredgenesignatureS1(1,:),ZscoredgenesignatureS2(1,:));
d=corrcoef(ZscoredgenesignatureS1(1,abs(ZscoredgenesignatureS1(1,:))>0.07),ZscoredgenesignatureS2(1,abs(ZscoredgenesignatureS1(1,:))>0.07));
title(strcat("well 1 vs well 2, corr1 ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,2)
scatter(ZscoredgenesignatureS1(1,:),ZscoredgenesignatureS3(1,:),'.')
c=corrcoef(ZscoredgenesignatureS1(1,:),ZscoredgenesignatureS3(1,:));
d=corrcoef(ZscoredgenesignatureS1(1,abs(ZscoredgenesignatureS1(1,:))>0.07),ZscoredgenesignatureS3(1,abs(ZscoredgenesignatureS1(1,:))>0.07));
title(strcat("well 1 vs well 3, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,3)
 scatter(ZscoredgenesignatureS1(1,:),ZscoredgenesignatureS4(1,:),'.')
c=corrcoef(ZscoredgenesignatureS1(1,:),ZscoredgenesignatureS4(1,:));
d=corrcoef(ZscoredgenesignatureS1(1,abs(ZscoredgenesignatureS1(1,:))>0.07),ZscoredgenesignatureS4(1,abs(ZscoredgenesignatureS1(1,:))>0.07));
title(strcat("well 1 vs well 4, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,4)
 scatter(ZscoredgenesignatureS1(1,:),ZscoredgenesignatureS5(1,:),'.')
c=corrcoef(ZscoredgenesignatureS1(1,:),ZscoredgenesignatureS5(1,:));
d=corrcoef(ZscoredgenesignatureS1(1,abs(ZscoredgenesignatureS1(1,:))>0.07),ZscoredgenesignatureS5(1,abs(ZscoredgenesignatureS1(1,:))>0.07));
title(strcat("well 1 vs well 5, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
sgtitle('zscored expression for cluster large peak, no norm') 
figure, 
subplot(2,2,1)
scatter(ZscoredgenesignatureS1(2,:),ZscoredgenesignatureS2(2,:),'.')
c=corrcoef(ZscoredgenesignatureS1(2,:),ZscoredgenesignatureS2(2,:));
d=corrcoef(ZscoredgenesignatureS1(2,abs(ZscoredgenesignatureS1(2,:))>0.07),ZscoredgenesignatureS2(2,abs(ZscoredgenesignatureS1(2,:))>0.07));
title(strcat("well 1 vs well 2, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,2)
scatter(ZscoredgenesignatureS1(2,:),ZscoredgenesignatureS3(2,:),'.')
c=corrcoef(ZscoredgenesignatureS1(2,:),ZscoredgenesignatureS3(2,:));
d=corrcoef(ZscoredgenesignatureS1(2,abs(ZscoredgenesignatureS1(2,:))>0.07),ZscoredgenesignatureS3(2,abs(ZscoredgenesignatureS1(2,:))>0.07));
title(strcat("well 1 vs well 3, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,3)
 scatter(ZscoredgenesignatureS1(2,:),ZscoredgenesignatureS4(2,:),'.')
c=corrcoef(ZscoredgenesignatureS1(2,:),ZscoredgenesignatureS4(2,:));
d=corrcoef(ZscoredgenesignatureS1(2,abs(ZscoredgenesignatureS1(2,:))>0.07),ZscoredgenesignatureS4(2,abs(ZscoredgenesignatureS1(2,:))>0.07));
title(strcat("well 1 vs well 4, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,4)
 scatter(ZscoredgenesignatureS1(2,:),ZscoredgenesignatureS5(2,:),'.')
c=corrcoef(ZscoredgenesignatureS1(2,:),ZscoredgenesignatureS5(2,:));
d=corrcoef(ZscoredgenesignatureS1(2,abs(ZscoredgenesignatureS1(2,:))>0.07),ZscoredgenesignatureS5(2,abs(ZscoredgenesignatureS1(2,:))>0.07));
title(strcat("well 1 vs well 5, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
sgtitle('zscored expression for cluster inactive cells, no norm') 
figure, 
subplot(2,2,1)
scatter(ZscoredgenesignatureS1(3,:),ZscoredgenesignatureS2(3,:),'.')
c=corrcoef(ZscoredgenesignatureS1(3,:),ZscoredgenesignatureS2(3,:));
d=corrcoef(ZscoredgenesignatureS1(3,abs(ZscoredgenesignatureS1(3,:))>0.07),ZscoredgenesignatureS2(3,abs(ZscoredgenesignatureS1(3,:))>0.07));
title(strcat("well 1 vs well 2, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,2)
scatter(ZscoredgenesignatureS1(3,:),ZscoredgenesignatureS3(3,:),'.')
c=corrcoef(ZscoredgenesignatureS1(3,:),ZscoredgenesignatureS3(3,:));
d=corrcoef(ZscoredgenesignatureS1(3,abs(ZscoredgenesignatureS1(3,:))>0.07),ZscoredgenesignatureS3(3,abs(ZscoredgenesignatureS1(3,:))>0.07));
title(strcat("well 1 vs well 3, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,3)
 scatter(ZscoredgenesignatureS1(3,:),ZscoredgenesignatureS4(3,:),'.')
c=corrcoef(ZscoredgenesignatureS1(3,:),ZscoredgenesignatureS4(3,:));
d=corrcoef(ZscoredgenesignatureS1(3,abs(ZscoredgenesignatureS1(3,:))>0.07),ZscoredgenesignatureS4(3,abs(ZscoredgenesignatureS1(3,:))>0.07));
title(strcat("well 1 vs well 4, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,4)
 scatter(ZscoredgenesignatureS1(3,:),ZscoredgenesignatureS5(3,:),'.')
c=corrcoef(ZscoredgenesignatureS1(3,:),ZscoredgenesignatureS5(3,:));
d=corrcoef(ZscoredgenesignatureS1(3,abs(ZscoredgenesignatureS1(3,:))>0.07),ZscoredgenesignatureS5(3,abs(ZscoredgenesignatureS1(3,:))>0.07));
title(strcat("well 1 vs well 5, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
sgtitle('zscored expression for cluster large early transients, no norm') 
figure, 
subplot(2,2,1)
scatter(ZscoredgenesignatureS1(4,:),ZscoredgenesignatureS2(4,:),'.')
c=corrcoef(ZscoredgenesignatureS1(4,:),ZscoredgenesignatureS2(4,:));
d=corrcoef(ZscoredgenesignatureS1(4,abs(ZscoredgenesignatureS1(4,:))>0.07),ZscoredgenesignatureS2(4,abs(ZscoredgenesignatureS1(4,:))>0.07));
title(strcat("well 1 vs well 2, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,2)
scatter(ZscoredgenesignatureS1(4,:),ZscoredgenesignatureS3(4,:),'.')
c=corrcoef(ZscoredgenesignatureS1(4,:),ZscoredgenesignatureS3(4,:));
d=corrcoef(ZscoredgenesignatureS1(4,abs(ZscoredgenesignatureS1(4,:))>0.07),ZscoredgenesignatureS3(4,abs(ZscoredgenesignatureS1(4,:))>0.07));
title(strcat("well 1 vs well 3, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,3)
 scatter(ZscoredgenesignatureS1(4,:),ZscoredgenesignatureS4(4,:),'.')
c=corrcoef(ZscoredgenesignatureS1(4,:),ZscoredgenesignatureS4(4,:));
d=corrcoef(ZscoredgenesignatureS1(4,abs(ZscoredgenesignatureS1(4,:))>0.07),ZscoredgenesignatureS4(4,abs(ZscoredgenesignatureS1(4,:))>0.07));
title(strcat("well 1 vs well 4, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,4)
 scatter(ZscoredgenesignatureS1(4,:),ZscoredgenesignatureS5(4,:),'.')
c=corrcoef(ZscoredgenesignatureS1(4,:),ZscoredgenesignatureS5(4,:));
d=corrcoef(ZscoredgenesignatureS1(4,abs(ZscoredgenesignatureS1(4,:))>0.07),ZscoredgenesignatureS2(4,abs(ZscoredgenesignatureS1(4,:))>0.07));
title(strcat("well 1 vs well 5, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
sgtitle('zscored expression for cluster small peak, no norm') 


figure, 
subplot(2,2,1)
scatter(ZscoredgenesignatureS1(5,:),ZscoredgenesignatureS2(5,:),'.')
c=corrcoef(ZscoredgenesignatureS1(5,:),ZscoredgenesignatureS2(5,:));
d=corrcoef(ZscoredgenesignatureS1(5,abs(ZscoredgenesignatureS1(5,:))>0.07),ZscoredgenesignatureS2(5,abs(ZscoredgenesignatureS1(5,:))>0.07));
title(strcat("well 1 vs well 2, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,2)
scatter(ZscoredgenesignatureS1(5,:),ZscoredgenesignatureS3(5,:),'.')
c=corrcoef(ZscoredgenesignatureS1(5,:),ZscoredgenesignatureS3(5,:));
d=corrcoef(ZscoredgenesignatureS1(5,abs(ZscoredgenesignatureS1(5,:))>0.07),ZscoredgenesignatureS3(5,abs(ZscoredgenesignatureS1(5,:))>0.07));
title(strcat("well 1 vs well 3, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,3)
 scatter(ZscoredgenesignatureS1(5,:),ZscoredgenesignatureS4(5,:),'.')
c=corrcoef(ZscoredgenesignatureS1(5,:),ZscoredgenesignatureS4(5,:));
d=corrcoef(ZscoredgenesignatureS1(5,abs(ZscoredgenesignatureS1(5,:))>0.07),ZscoredgenesignatureS4(5,abs(ZscoredgenesignatureS1(5,:))>0.07));
title(strcat("well 1 vs well 4, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,4)
 scatter(ZscoredgenesignatureS1(5,:),ZscoredgenesignatureS5(5,:),'.')
c=corrcoef(ZscoredgenesignatureS1(5,:),ZscoredgenesignatureS5(5,:));
d=corrcoef(ZscoredgenesignatureS1(5,abs(ZscoredgenesignatureS1(5,:))>0.07),ZscoredgenesignatureS5(5,abs(ZscoredgenesignatureS1(5,:))>0.07));
title(strcat("well 1 vs well 5, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
sgtitle('zscored expression for cluster step, no norm') 

figure, 
subplot(2,2,1)
scatter(ZscoredgenesignatureS1(6,:),ZscoredgenesignatureS2(6,:),'.')
c=corrcoef(ZscoredgenesignatureS1(6,:),ZscoredgenesignatureS2(6,:));
d=corrcoef(ZscoredgenesignatureS1(6,abs(ZscoredgenesignatureS1(6,:))>0.07),ZscoredgenesignatureS2(6,abs(ZscoredgenesignatureS1(6,:))>0.07));
title(strcat("well 1 vs well 2, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,2)
scatter(ZscoredgenesignatureS1(6,:),ZscoredgenesignatureS3(6,:),'.')
c=corrcoef(ZscoredgenesignatureS1(6,:),ZscoredgenesignatureS3(6,:));
d=corrcoef(ZscoredgenesignatureS1(6,abs(ZscoredgenesignatureS1(6,:))>0.07),ZscoredgenesignatureS3(6,abs(ZscoredgenesignatureS1(6,:))>0.07));
title(strcat("well 1 vs well 3, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,3)
 scatter(ZscoredgenesignatureS1(6,:),ZscoredgenesignatureS4(6,:),'.')
c=corrcoef(ZscoredgenesignatureS1(6,:),ZscoredgenesignatureS4(6,:));
d=corrcoef(ZscoredgenesignatureS1(6,abs(ZscoredgenesignatureS1(6,:))>0.07),ZscoredgenesignatureS4(6,abs(ZscoredgenesignatureS1(6,:))>0.07));
title(strcat("well 1 vs well 4, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
subplot(2,2,4)
 scatter(ZscoredgenesignatureS1(6,:),ZscoredgenesignatureS5(6,:),'.')
c=corrcoef(ZscoredgenesignatureS1(6,:),ZscoredgenesignatureS5(6,:));
d=corrcoef(ZscoredgenesignatureS1(6,abs(ZscoredgenesignatureS1(6,:))>0.07),ZscoredgenesignatureS5(6,abs(ZscoredgenesignatureS1(6,:))>0.07));
title(strcat("well 1 vs well 5, correlation ",num2str(c(1,2))," corr2 ",num2str(d(1,2))))
sgtitle('zscored expression for cluster delayed, no norm')
%% see the traces of the 6 clusters
allcalciumtracesnormed=table2array(readtable('\\helium\broad_clearylab\Users\Loic\pooledrevision\allcalciumtracesnormed_0329.csv'));
idx=table2array(readtable('\\helium\broad_clearylab\Users\Loic\pooledrevision\singleclusters\calciumclusterall_idx_it_calciumZscoreONEatAtime.csv'));
idx(:,4)=idx(:,3);
idx(idx(:,4)==7,4)=5;
idx(idx(:,4)==6,4)=2;
idx(idx(:,4)==8,4)=6;
clusternames=["large peak","inactive","large early transient","small peak","step","delayed"];

calcium=allcalciumtracesnormed;
time=[1:200];
for sample=1:5
        if sample==1
            thissampleID=21;
        elseif sample==2
            thissampleID=22;
        elseif sample==3
            thissampleID=31;          
        elseif sample==4
            thissampleID=32;           
        elseif sample==5
            thissampleID=4;
        end
        figure,
for i=1:6
    subtraces=calcium(idx(:,3)==i & idx(:,1)==thissampleID,:);
    for j=1:10
        subplot(2,3,i)
        plot(time,subtraces(10+j,:)+(j)*20,'LineWidth',2)
    %     title(clusternames(i))
        ylim([0 350])
        axis off
        hold on
    end
    title(strcat("Cluster ",clusternames(i)))
    sgtitle(strcat('Traces for sample ',num2str(sample)))

end
set(gcf,'color','w');
end
