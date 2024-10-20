uiopen('\\helium\broad_clearylab\Users\Loic\FRonTumors\tumordata\split1_LFCstake1.csv',1)
uiopen('\\helium\broad_clearylab\Users\Loic\FRonTumors\tumordata\split1_qvalstake1.csv',1)
uiopen('\\helium\broad_clearylab\Users\Loic\FRonTumors\tumordata\split2_LFCstake1.csv',1)
uiopen('\\helium\broad_clearylab\Users\Loic\FRonTumors\tumordata\split2_qvalstake1.csv',1)

lfc1=split1LFCstake1(:,1:35)';
qval1=split1qvalstake1(:,1:35)';
lfc2=split2LFCstake1(:,1:35)';
qval2=split2qvalstake1(:,1:35)';

oldlfc1=lfc1;
oldlfc2=lfc2;
oldqval1=qval1;
oldqval2=qval2;

perturbationNames=perturbationNames(1:35);
flfc1=reshape(lfc1,[size(lfc1,1)*size(lfc1,2),1]);
fqval1=reshape(qval1,[size(lfc1,1)*size(lfc1,2),1]);
flfc2=reshape(lfc2,[size(lfc1,1)*size(lfc1,2),1]);
fqval2=reshape(qval2,[size(lfc1,1)*size(lfc1,2),1]);
keep=[];
for i=1:size(lfc1,2)
    keep=[keep;keepguides];
end
figure, scatter(flfc1(keep>0&(fqval1<0.1|fqval2<0.1)),flfc2(keep>0&(fqval1<0.1|fqval2<0.1)),'.')
hold on
scatter(flfc1(keep>0&(fqval1<0.1&fqval2<0.1)),flfc2(keep>0&(fqval1<0.1&fqval2<0.1)),'.')

for i=1:35
    figure,
    X=lfc1(i,qval1(i,:)<0.1&qval2(i,:)<0.1);
    Y=lfc2(i,qval1(i,:)<0.1&qval2(i,:)<0.1);
    scatter(X,Y,'.')
    xlim([-1,1])
    ylim([-1,1])
    title(perturbationNames(i))
end
xlabel("Split 1")
ylabel("Split 2")
set(gcf,'color','w');caxis([-1 1])
correlationeither=corrcoef(flfc1(keep>0&(fqval1<0.1|fqval2<0.1)),flfc2(keep>0&(fqval1<0.1|fqval2<0.1)));
correlationboth=corrcoef(flfc1(keep>0&(fqval1<0.1&fqval2<0.1)),flfc2(keep>0&(fqval1<0.1&fqval2<0.1)));
title(strcat("Significant effects from 2 random splits, correlations:", num2str(correlationeither(1,2)), " and " , num2str(correlationboth(1,2))))

figure, heatmap(geneNames,perturbationNames,lfc1,'Colormap',redblue(256))
set(gcf,'color','w');caxis([-0.8 0.8])
figure, heatmap(geneNames,perturbationNames,lfc2,'Colormap',redblue(256))
set(gcf,'color','w');caxis([-0.8 0.8])


%% check between guides
uiopen('\\helium\broad_clearylab\Users\Loic\FRonTumors\tumordata\notpoolednot2_LFCs.csv',1)
uiopen('\\helium\broad_clearylab\Users\Loic\FRonTumors\tumordata\normalno2notpooled_qvals.csv',1)
guides=perturbationNames;
for i=1:35
    figure
    scatter(notpoolednot2LFCs(notpoolednot2qvals(:,2*i-1)<0.1|notpoolednot2qvals(:,2*i)<0.1,2*i-1),notpoolednot2LFCs(notpoolednot2qvals(:,2*i-1)<0.1|notpoolednot2qvals(:,2*i)<0.1,2*i),'.')
    hold on
    scatter(notpoolednot2LFCs(notpoolednot2qvals(:,2*i-1)<0.1&notpoolednot2qvals(:,2*i)<0.1,2*i-1),notpoolednot2LFCs(notpoolednot2qvals(:,2*i-1)<0.1&notpoolednot2qvals(:,2*i)<0.1,2*i),'.')
    c=corrcoef(notpoolednot2LFCs(notpoolednot2qvals(:,2*i-1)<0.1& notpoolednot2qvals(:,2*i)<0.1,2*i-1),notpoolednot2LFCs(notpoolednot2qvals(:,2*i-1)<0.1&notpoolednot2qvals(:,2*i)<0.1,2*i)); 
    if size(c,2)>1&~isnan(c(1,1))
        title(strcat(guides(2*i), ' correlation', num2str(c(1,2))))
        correlations(i)=c(1,2);
    else
%         title(strcat(guides(2*i), ' correlation', num2str(c(1,1))))
%         correlations(i)=c(1,1);
    end
end
for i=1:35
    guides1(:,i)=notpoolednot2LFCs(:,2*i-1);
    guides2(:,i)=notpoolednot2LFCs(:,2*i);
end
writematrix(guides1',fullfile(mypath,'guides1effects.csv'))
writematrix(guides2',fullfile(mypath,'guides2effects.csv'))



%%
load('\\helium\broad_clearylab\Users\Loic\nicetables\workspace0903_1400.mat')
filteredsplit1LFCs=lfc1(keepguides,keepgenes);
filteredsplit2LFCs=lfc2(keepguides,keepgenes);
figure, heatmap(geneNames(keepgenes),perturbationNames(keepguides),filteredsplit1LFCs,'Colormap',redblue(256))
set(gcf,'color','w');caxis([-0.5 0.5])
whitened=filterednormalLFCs;whitened(filterednormalqvals>0.1)=0;
figure, heatmap(geneNames(keepgenes),perturbationNames(keepguides),filteredsplit2LFCs,'Colormap',redblue(256))
set(gcf,'color','w');caxis([-0.5 0.5])


tosortsplit1=filteredsplit1LFCs;
sortedSplit1=tosortsplit1(idxLFCguides==1,:);
for i=2:max(idxLFCguides)
sortedSplit1=[sortedSplit1;tosortsplit1(idxLFCguides==i,:)];
end
tosortsplit1=sortedSplit1';
% [idxLFCgenes,C] = kmeans(tosort,5,'Distance','correlation','Replicates',60,'MaxIter',1000,'Options',opts);
sortedSplit1=tosortsplit1(idxLFCgenes==1,:);
for i=2:max(idxLFCgenes)
sortedSplit1=[sortedSplit1;tosortsplit1(idxLFCgenes==i,:)];
end

tosortsplit2=filteredsplit2LFCs;
sortedSplit2=tosortsplit2(idxLFCguides==1,:);
for i=2:max(idxLFCguides)
sortedSplit2=[sortedSplit2;tosortsplit2(idxLFCguides==i,:)];
end
tosortsplit2=sortedSplit2';
% [idxLFCgenes,C] = kmeans(tosort,5,'Distance','correlation','Replicates',60,'MaxIter',1000,'Options',opts);
sortedSplit2=tosortsplit1(idxLFCgenes==1,:);
for i=2:max(idxLFCgenes)
sortedSplit2=[sortedSplit2;tosortsplit2(idxLFCgenes==i,:)];
end
writematrix(sortedSplit1',fullfile(mypath,'sortedsplit1.csv'))
writematrix(sortedSplit2',fullfile(mypath,'sortedsplit2.csv'))

figure, heatmap(sortedgenes,sortedguides,sortedSplit1','Colormap',redblue(256))
xlabel("Gene")
ylabel("Perturbation")
set(gcf,'color','w');caxis([-0.7 0.7])
title("Significant effects without immune neighbor")

figure, heatmap(sortedgenes,sortedguides,sortedSplit2','Colormap',redblue(256))
xlabel("Gene")
ylabel("Perturbation")
set(gcf,'color','w');caxis([-0.7 0.7])
title("Significant effects without immune neighbor")
%%
lfc1=oldlfc1;
lfc2=oldlfc2;
qval1=oldqval1;
qval2=oldqval2;
filteredsplit1LFCs=lfc1(keepguides,keepgenes);
filteredsplit2LFCs=lfc2(keepguides,keepgenes);

tosortsplit1=filteredsplit1LFCs;
sortedSplit1=tosortsplit1(idxLFCguides==1,:);
for i=2:max(idxLFCguides)
sortedSplit1=[sortedSplit1;tosortsplit1(idxLFCguides==i,:)];
end
tosortsplit1=sortedSplit1';
% [idxLFCgenes,C] = kmeans(tosort,5,'Distance','correlation','Replicates',60,'MaxIter',1000,'Options',opts);
sortedSplit1=tosortsplit1(idxLFCgenes==1,:);
for i=2:max(idxLFCgenes)
sortedSplit1=[sortedSplit1;tosortsplit1(idxLFCgenes==i,:)];
end

tosortsplit2=filteredsplit2LFCs;
sortedSplit2=tosortsplit2(idxLFCguides==1,:);
for i=2:max(idxLFCguides)
sortedSplit2=[sortedSplit2;tosortsplit2(idxLFCguides==i,:)];
end
tosortsplit2=sortedSplit2';
% [idxLFCgenes,C] = kmeans(tosort,5,'Distance','correlation','Replicates',60,'MaxIter',1000,'Options',opts);
sortedSplit2=tosortsplit1(idxLFCgenes==1,:);
for i=2:max(idxLFCgenes)
sortedSplit2=[sortedSplit2;tosortsplit2(idxLFCgenes==i,:)];
end
writematrix(sortedSplit1',fullfile(mypath,'sortedsplit1.csv'))
writematrix(sortedSplit2',fullfile(mypath,'sortedsplit2.csv'))

figure, heatmap(sortedgenes,sortedguides,sortedSplit1','Colormap',redblue(256))
xlabel("Gene")
ylabel("Perturbation")
set(gcf,'color','w');caxis([-0.7 0.7])
title("Significant effects without immune neighbor")

figure, heatmap(sortedgenes,sortedguides,sortedSplit2','Colormap',redblue(256))
xlabel("Gene")
ylabel("Perturbation")
set(gcf,'color','w');caxis([-0.7 0.7])
title("Significant effects without immune neighbor")

