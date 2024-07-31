perturbcalcium=table2array(readtable('\\helium\broad_clearylab\Users\Loic\AsdRevstables\clusteredCalciumperturbcells.csv'));
perturbations=table2array(readtable('\\helium\broad_clearylab\Users\Loic\AsdRevstables\perturbationsInCalcium.csv'));

perturbcalciumlessclusters=perturbcalcium;
perturbcalciumlessclusters(perturbcalciumlessclusters(:,3)==7,3)=5;
perturbcalciumlessclusters(perturbcalciumlessclusters(:,3)==6,3)=2;
perturbcalciumlessclusters(perturbcalciumlessclusters(:,3)==8,3)=6;

guidespercluster=zeros(6,128);
for cluster=1:max(perturbcalciumlessclusters(:,3))
    theseperturbations=perturbations(perturbcalciumlessclusters(:,3)==cluster,:);
    mylist=theseperturbations(find(theseperturbations));
    for i=1:128
        guidespercluster(cluster,i)=sum(mylist==i);
    end
end
cellspercluster=[]
for i=1:max(perturbcalciumlessclusters(:,3))
cellspercluster=[cellspercluster;sum(perturbcalciumlessclusters(:,3)==i)];
end
expectedfreq=zeros(6,128);
for j=1:128
    expectedfreq(:,j)=sum(guidespercluster(:,j))*cellspercluster'/sum(cellspercluster);
end
FC=log(guidespercluster./expectedfreq);
clusternames=["large peak","inactive","large early transient","small peak","step","delayed"];
% heatmap(guides(1:128),clusternames,FC,'Colormap',redblue(256))
% ylabel("Cluster")
% xlabel("Perturbation")
% title("Perturbation enrichment per cluster")
set(gcf,'color','w');caxis([-1.5 1.5])

M=size(perturbations,1);%sample size
K=sum(guidespercluster(1,:));%number of good items in pop
N=sum(guidespercluster(:,1));%number drawn
pval=1-hygecdf(guidespercluster(:,1)-1,M,K,N);
p=zeros(max(perturbcalciumlessclusters(:,3)),size(guidespercluster,2));
for i=1:size(guidespercluster,2)
    for j=1:size(guidespercluster,1)
        K=sum(guidespercluster(j,:));
        N=sum(guidespercluster(:,i));
        pval=min(1-hygecdf(guidespercluster(j,i)-1,M,K,N),hygecdf(guidespercluster(j,i)-1,M,K,N));
        p(j,i)=pval;
    end
end
globalFC=FC;
globalp=p;
% sigEnriched=sum(globalp<0.1,1);
% sigEnriched=[sigEnriched;sigEnriched;sigEnriched;sigEnriched;sigEnriched;sigEnriched];
%% 
figure,
    idx = randperm(size(perturbations,1)) ;
shuffled = perturbations(idx,:) ;
shuffledcluster=perturbcalciumlessclusters(idx,:) ;
split1=shuffled(1:floor(size(shuffled,1)/2),:);
split2=shuffled(floor(size(shuffled,1)/2):end,:);
split1cluster=shuffledcluster(1:floor(size(shuffled,1)/2),:);
split2cluster=shuffledcluster(floor(size(shuffled,1)/2):end,:);

%% split 1
guidespercluster=zeros(6,128);
for cluster=1:max(split1cluster(:,3))
    theseperturbations=split1(split1cluster(:,3)==cluster,:);
    mylist=theseperturbations(find(theseperturbations));
    for i=1:128
        guidespercluster(cluster,i)=sum(mylist==i);
    end
end

cellspercluster=[]
for i=1:max(split1cluster(:,3))
cellspercluster=[cellspercluster;sum(split1cluster(:,3)==i)];
end
expectedfreq=zeros(6,128);
for j=1:128
    expectedfreq(:,j)=sum(guidespercluster(:,j))*cellspercluster'/sum(cellspercluster);
end
FC=log(guidespercluster./expectedfreq);
clusternames=["large peak","inactive","large early transient","small peak","step","delayed"];
% figure, heatmap(guides(1:128),clusternames,FC,'Colormap',redblue(256))
% ylabel("Cluster")
% xlabel("Perturbation")
% title("Perturbation enrichment per cluster, split2")
% set(gcf,'color','w');caxis([-1.5 1.5])
FCsplit1=FC;

M=size(split1,1);%sample size
K=sum(guidespercluster(1,:));%number of good items in pop
N=sum(guidespercluster(:,1));%number drawn
pval=1-hygecdf(guidespercluster(:,1)-1,M,K,N);
p=zeros(max(perturbcalciumlessclusters(:,3)),size(guidespercluster,2));
for i=1:size(guidespercluster,2)
    for j=1:size(guidespercluster,1)
        K=sum(guidespercluster(j,:));
        N=sum(guidespercluster(:,i));
        pval=min(1-hygecdf(guidespercluster(j,i)-1,M,K,N),hygecdf(guidespercluster(j,i)-1,M,K,N));
        p(j,i)=pval;
    end
end
psplit1=p;
%% split 2
guidespercluster=zeros(6,128);
for cluster=1:max(split2cluster(:,3))
    theseperturbations=split2(split2cluster(:,3)==cluster,:);
    mylist=theseperturbations(find(theseperturbations));
    for i=1:128
        guidespercluster(cluster,i)=sum(mylist==i);
    end
end

cellspercluster=[]
for i=1:max(split2cluster(:,3))
cellspercluster=[cellspercluster;sum(split2cluster(:,3)==i)];
end
expectedfreq=zeros(6,128);
for j=1:128
    expectedfreq(:,j)=sum(guidespercluster(:,j))*cellspercluster'/sum(cellspercluster);
end
FC=log(guidespercluster./expectedfreq);

M=size(split2,1);%sample size
K=sum(guidespercluster(1,:));%number of good items in pop
N=sum(guidespercluster(:,1));%number drawn
pval=1-hygecdf(guidespercluster(:,1)-1,M,K,N);
p=zeros(max(perturbcalciumlessclusters(:,3)),size(guidespercluster,2));
for i=1:size(guidespercluster,2)
    for j=1:size(guidespercluster,1)
        K=sum(guidespercluster(j,:));
        N=sum(guidespercluster(:,i));
        pval=min(1-hygecdf(guidespercluster(j,i)-1,M,K,N),hygecdf(guidespercluster(j,i)-1,M,K,N));
        p(j,i)=pval;
    end
end
psplit2=p;
clusternames=["large peak","inactive","large early transient","small peak","step","delayed"];
% figure, heatmap(guides(1:128),clusternames,FC,'Colormap',redblue(256))
% ylabel("Cluster")
% xlabel("Perturbation")
% title("Perturbation enrichment per cluster, split2")
% set(gcf,'color','w');caxis([-1.5 1.5])
flatp1=reshape(psplit1,[1, size(FCsplit1,1)*size(FCsplit1,2)]);
flatp2=reshape(psplit2,[1, size(FC,1)*size(FC,2)]);
flatFCsplit1=reshape(FCsplit1,[1, size(FCsplit1,1)*size(FCsplit1,2)]);
flatFCsplit2=reshape(FC,[1, size(FC,1)*size(FC,2)]);
flatglobalFC=reshape(globalFC,[1, size(FC,1)*size(FC,2)]);;
flatglobalp=reshape(globalp,[1, size(FC,1)*size(FC,2)]);;

% subplot(2,2,2)
% 
% scatter(flatFCsplit1(psplit1<0.05&psplit2<0.05),flatFCsplit2(psplit1<0.05&psplit2<0.05),'.')
% 
% scatter(flatFCsplit1(globalp<0.1),globalFC(globalp<0.1),'.')
% hold on
% scatter(flatFCsplit2(globalp<0.1),globalFC(globalp<0.1),'.')
% hold on
flatFCsplit1(isinf(flatFCsplit1))=sign(flatFCsplit1(isinf(flatFCsplit1)))*2;
flatFCsplit2(isinf(flatFCsplit2))=sign(flatFCsplit2(isinf(flatFCsplit2)))*2;
flatglobalFC(isinf(flatglobalFC))=sign(flatglobalFC(isinf(flatglobalFC)))*2;
subplot(2,2,1)
scatter(flatFCsplit1(flatp1<0.05&flatp2<0.05 ),flatFCsplit2(flatp1<0.05&flatp2<0.05),'.')
title("Guide enrichment/depletion significant in both halves")
xlabel('Split 1')
ylabel('Split 2')
C=corrcoef(flatFCsplit1(flatp1<0.05&flatp2<0.05&~isinf(flatFCsplit1)&~isinf(flatFCsplit2)),flatFCsplit2(flatp1<0.05&flatp2<0.05&~isinf(flatFCsplit1)&~isinf(flatFCsplit2)));
annotation('textbox',[0.15 0.9 0 0],'string',strcat('r=',num2str(C(1,2))),'FitBoxToText','on','EdgeColor','black')     

subplot(2,2,2)
scatter(flatFCsplit1(psplit1<0.05 | psplit2<0.05),flatFCsplit2(psplit1<0.05 | psplit2<0.05),'.')
title("Guide enrichment/depletion significant in at least one half")
xlabel('Split 1')
ylabel('Split 2')
C=corrcoef(flatFCsplit1(psplit1<0.05 | psplit2<0.05),flatFCsplit2(psplit1<0.05 | psplit2<0.05));
annotation('textbox',[0.6 0.9 0 0],'string',strcat('r=',num2str(C(1,2))),'FitBoxToText','on','EdgeColor','black')     

subplot(2,2,3)
scatter(flatglobalFC(flatglobalp<0.05),flatFCsplit1(flatglobalp<0.05),'.')
title("Guide enrichment/depletion significant in the pooled data")
xlabel('Pooled data')
ylabel('Split 1')
C=corrcoef(flatFCsplit1(flatglobalp<0.05),flatglobalFC(flatglobalp<0.05));
annotation('textbox',[0.15 0.45 0 0],'string',strcat('r=',num2str(C(1,2))),'FitBoxToText','on','EdgeColor','black')     


subplot(2,2,4)
scatter(flatglobalFC(flatglobalp<0.05),flatFCsplit2(flatglobalp<0.05),'.')
title("Guide enrichment/depletion significant in the pooled data")
xlabel('Pooled data')
ylabel('Split 2')
set(gcf,'color','w');
C=corrcoef(flatFCsplit2(flatglobalp<0.05),flatglobalFC(flatglobalp<0.05));
annotation('textbox',[0.6 0.45 0 0],'string',strcat('r=',num2str(C(1,2))),'FitBoxToText','on','EdgeColor','black')     

