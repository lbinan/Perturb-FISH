calciumcluster=table2array(readtable('\\helium\broad_clearylab\Users\Loic\pooledrevision\singleclusters\calciumclusterall_idx_it_calciumZscoreONEatAtime.csv'));
pooledmerfish=table2array(readtable('\\helium\broad_clearylab\Users\Loic\AsdRevstables\pooledmerfishmatchedtocalcium.csv'));
%load all genes
%makes heatmap of gene zscored signature per calcium phenotype, with genes
%sorted in the same order as the clustering from FR infered perturbation
%effects
pooledmerfish=pooledmerfish(:,3:end);
pooledmerfish=pooledmerfish(:,1:485);
totalcounts=sum(pooledmerfish,2);
badcells=totalcounts<200;
badcells=badcells|(totalcounts>1900);
pooledmerfish=pooledmerfish(~badcells,:);
for i=1:size(genes,1)
    if sum(strcmp(allgenes,genes(i)))>0
    merfish(:,i)=pooledmerfish(:,strcmp(allgenes,genes(i)));
    end
end
for i=1:size(keptgenes,2)
    if sum(strcmp(allgenes,keptgenes(i)))>0
    merfishkeptgenes(:,i)=pooledmerfish(:,strcmp(allgenes,keptgenes(i)));
    end
end
calcium=calciumcluster(~badcells,:);
zscored=merfishkeptgenes;
zscored(calcium(:,1)==21,:)=(merfishkeptgenes(calcium(:,1)==21,:)-mean(merfishkeptgenes(calcium(:,1)==21,:),1))./std(merfishkeptgenes(calcium(:,1)==21,:),1);
zscored(calcium(:,1)==22,:)=(merfishkeptgenes(calcium(:,1)==22,:)-mean(merfishkeptgenes(calcium(:,1)==22,:),1))./std(merfishkeptgenes(calcium(:,1)==22,:),1);
zscored(calcium(:,1)==31,:)=(merfishkeptgenes(calcium(:,1)==31,:)-mean(merfishkeptgenes(calcium(:,1)==31,:),1))./std(merfishkeptgenes(calcium(:,1)==31,:),1);
zscored(calcium(:,1)==32,:)=(merfishkeptgenes(calcium(:,1)==32,:)-mean(merfishkeptgenes(calcium(:,1)==32,:),1))./std(merfishkeptgenes(calcium(:,1)==32,:),1);
zscored(calcium(:,1)==4,:)=(merfishkeptgenes(calcium(:,1)==4,:)-mean(merfishkeptgenes(calcium(:,1)==4,:),1))./std(merfishkeptgenes(calcium(:,1)==4,:),1);


calciumlessclusters=calcium;
calciumlessclusters(calciumlessclusters(:,3)==7,3)=5;
calciumlessclusters(calciumlessclusters(:,3)==6,3)=2;
calciumlessclusters(calciumlessclusters(:,3)==8,3)=6;

for j=1:max(calciumlessclusters(:,3))
    Zscoredgenesignature(j,:)=mean(zscored(calciumlessclusters(:,3)==j,:));
end
sortedZscore=Zscoredgenesignature(:,idxgenes==1);
for i=2:max(idxgenes)
sortedZscore=[sortedZscore,Zscoredgenesignature(:,idxgenes==i)];
end
figure, heatmap(sortedgenes,clusternames,sortedZscore,'Colormap',redblue(256))
set(gcf,'color','w');
caxis([-0.3 0.3])
%% normed to total counts
normedmerfish=merfishkeptgenes./sum(merfishkeptgenes,2);
normedzscored=normedmerfish;
normedzscored(calcium(:,1)==21,:)=(normedmerfish(calcium(:,1)==21,:)-mean(normedmerfish(calcium(:,1)==21,:),1))./std(normedmerfish(calcium(:,1)==21,:),1);
normedzscored(calcium(:,1)==22,:)=(normedmerfish(calcium(:,1)==22,:)-mean(normedmerfish(calcium(:,1)==22,:),1))./std(normedmerfish(calcium(:,1)==22,:),1);
normedzscored(calcium(:,1)==31,:)=(normedmerfish(calcium(:,1)==31,:)-mean(normedmerfish(calcium(:,1)==31,:),1))./std(normedmerfish(calcium(:,1)==31,:),1);
normedzscored(calcium(:,1)==32,:)=(normedmerfish(calcium(:,1)==32,:)-mean(normedmerfish(calcium(:,1)==32,:),1))./std(normedmerfish(calcium(:,1)==32,:),1);
normedzscored(calcium(:,1)==4,:)=(normedmerfish(calcium(:,1)==4,:)-mean(normedmerfish(calcium(:,1)==4,:),1))./std(normedmerfish(calcium(:,1)==4,:),1);
for j=1:max(calciumlessclusters(:,3))
    normedZscoredgenesignature(j,:)=mean(normedzscored(calciumlessclusters(:,3)==j,:));
end
normedsortedZscore=normedZscoredgenesignature(:,idxgenes==1);
for i=2:max(idxgenes)
normedsortedZscore=[normedsortedZscore,normedZscoredgenesignature(:,idxgenes==i)];
end
figure, heatmap(sortedgenes,clusternames,normedsortedZscore,'Colormap',redblue(256))
set(gcf,'color','w');
caxis([-0.3 0.3])
for i=1:max(idxgenes)
    genescore(:,i)=sum(merfishkeptgenes(:,idxgenes==i),2);
end
genescore=genescore;
genescore(calcium(:,1)==21,:)=(genescore(calcium(:,1)==21,:)-mean(genescore(calcium(:,1)==21,:),1))./std(genescore(calcium(:,1)==21,:),1);
genescore(calcium(:,1)==22,:)=(genescore(calcium(:,1)==22,:)-mean(genescore(calcium(:,1)==22,:),1))./std(genescore(calcium(:,1)==22,:),1);
genescore(calcium(:,1)==31,:)=(genescore(calcium(:,1)==31,:)-mean(genescore(calcium(:,1)==31,:),1))./std(genescore(calcium(:,1)==31,:),1);
genescore(calcium(:,1)==32,:)=(genescore(calcium(:,1)==32,:)-mean(genescore(calcium(:,1)==32,:),1))./std(genescore(calcium(:,1)==32,:),1);
genescore(calcium(:,1)==4,:)=(genescore(calcium(:,1)==4,:)-mean(genescore(calcium(:,1)==4,:),1))./std(genescore(calcium(:,1)==4,:),1);
for j=1:max(calciumlessclusters(:,3))
    clustersignature(j,:)=mean(genescore(calciumlessclusters(:,3)==j,:));
end
figure, heatmap(["1","2","3","4","5"],clusternames,clustersignature,'Colormap',redblue(256),'CellLabelColor','none')
set(gcf,'color','w');
title("Z-scored expression of clutered genes")

for i=1:max(idxgenes)
    genescore(:,i)=sum(normedmerfish(:,idxgenes==i),2);
end
genescore=genescore;
genescore(calcium(:,1)==21,:)=(genescore(calcium(:,1)==21,:)-mean(genescore(calcium(:,1)==21,:),1))./std(genescore(calcium(:,1)==21,:),1);
genescore(calcium(:,1)==22,:)=(genescore(calcium(:,1)==22,:)-mean(genescore(calcium(:,1)==22,:),1))./std(genescore(calcium(:,1)==22,:),1);
genescore(calcium(:,1)==31,:)=(genescore(calcium(:,1)==31,:)-mean(genescore(calcium(:,1)==31,:),1))./std(genescore(calcium(:,1)==31,:),1);
genescore(calcium(:,1)==32,:)=(genescore(calcium(:,1)==32,:)-mean(genescore(calcium(:,1)==32,:),1))./std(genescore(calcium(:,1)==32,:),1);
genescore(calcium(:,1)==4,:)=(genescore(calcium(:,1)==4,:)-mean(genescore(calcium(:,1)==4,:),1))./std(genescore(calcium(:,1)==4,:),1);
for j=1:max(calciumlessclusters(:,3))
    clustersignature(j,:)=mean(genescore(calciumlessclusters(:,3)==j,:));
end
figure, heatmap(["1","2","3","4","5"],clusternames,clustersignature,'Colormap',redblue(256),'CellLabelColor','none')
set(gcf,'color','w');
title("Z-scored expression of clutered genes")

for j=1:max(idxgenes(:,1))
    meanclustersignature(:,j)=median(normedZscoredgenesignature(:,idxgenes(:,1)==j),2);
end
writematrix(normedsortedZscore,'\\helium\broad_clearylab\Users\Loic\AsdRevstables\calciumtomerfishMatching\sortedexpressionSignatures_normedtotalcounst.csv')
writematrix(sortedgenes,'\\helium\broad_clearylab\Users\Loic\AsdRevstables\calciumtomerfishMatching\sortedgenes.csv')

for i=1:142
    p(i)=anova1(normedmerfish(:,i),calcium(:,3),'off');
end
sortedp=p(:,idxgenes==1);
for i=2:max(idxgenes)
sortedp=[sortedp,p(:,idxgenes==i)];
end
writematrix(sortedp,'\\helium\broad_clearylab\Users\Loic\AsdRevstables\calciumtomerfishMatching\sortedANNOVAp.csv')


