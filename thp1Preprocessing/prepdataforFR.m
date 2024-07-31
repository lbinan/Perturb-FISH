%filter count tables, and reshape them to be the input that FR-perturb
%needs
lowcounts=sum(merfishtableonmosaicbacktomax(:,8:end),2)<70;%remove "cells" with low couts (i.e. debris)
filteredzombie=zombietableonmosaicbacktomax(~lowcounts,:);
filteredmerfih=merfishtableonmosaicbacktomax(~lowcounts,:);

highcounts=sum(filteredmerfih(:,8:end),2)>2300;%remove cells with high counts (i.e. doublets)
filteredzombie=filteredzombie(~highcounts,:);
filteredmerfih=filteredmerfih(~highcounts,:);


toobig=filteredmerfih(:,3)>80000;%remove cells that are too big (i.e. doublets)

filteredagainzombie=filteredzombie(~toobig,:);
filteredagainmerfish=filteredmerfih(~toobig,:);
% toosmall=filteredagainmerfish(:,3)<6000;

% filteredagainzombie=filteredagainzombie(~toosmall,:);
% filteredagainmerfish=filteredagainmerfish(~toosmall,:);
% 

justoneguide=filteredagainzombie(:,8)==1;%|filteredzombie(:,8)==2;%keep cells with only one guide
finalzombie=filteredagainzombie(justoneguide,:);
finalmerfish=filteredagainmerfish(justoneguide,:);
lpsminusindeces=ismember(finalzombie(:,1),listoflpsminusfovs);%remov cells the LPS minus well
finalzombie=finalzombie(~lpsminusindeces,:);
finalmerfish=finalmerfish(~lpsminusindeces,:);
%% keep only controls that only have one guides
onlyoneguideindeces=finalzombie(:,8)==1;
controlguides=sum(finalzombie(:,(8+71):(8+74)),2)>0;
finalmerfish=finalmerfish(~controlguides|(controlguides & onlyoneguideindeces),:);
finalzombie=finalzombie(~controlguides|(controlguides & onlyoneguideindeces),:);

%% save
covariates=finalmerfish(:,1:7);
merfishcounttable=finalmerfish(:,9:138);
zombiecountable=finalzombie(:,9:82);
writematrix(covariates','\\helium\broad_clearylab\Users\Loic\thp1homedish2_zombie\DataProcessedOnMosaic\covariatesTrTER.csv')
writematrix(covariates,'\\helium\broad_clearylab\Users\Loic\thp1homedish2_zombie\DataProcessedOnMosaic\covariatesTER.csv')
writematrix(merfishcounttable','\\helium\broad_clearylab\Users\Loic\thp1homedish2_zombie\DataProcessedOnMosaic\merfishcounttableTrTER.csv')
writematrix(merfishcounttable,'\\helium\broad_clearylab\Users\Loic\thp1homedish2_zombie\DataProcessedOnMosaic\merfishcounttableTER.csv')
writematrix(zombiecountable','\\helium\broad_clearylab\Users\Loic\thp1homedish2_zombie\DataProcessedOnMosaic\zombiecountableTrTER.csv')
writematrix(genelist','\\helium\broad_clearylab\Users\Loic\thp1homedish2_zombie\DataProcessedOnMosaic\genelistTER.csv')
writematrix(guidelist','\\helium\broad_clearylab\Users\Loic\thp1homedish2_zombie\DataProcessedOnMosaic\guidelistTER.csv')

writematrix(merfishcounttable,'\\helium\broad_clearylab\Users\Loic\thp1homedish2_zombie\DataProcessedOnMosaic\merfishcounttableTER.csv')

writematrix(covariates,'\\helium\broad_clearylab\Users\Loic\thp1homedish2_zombie\DataProcessedOnMosaic\covariatesTER.csv')
zombiecountableTr=zombiecountable';
for i=1:35
    pooledguides(i,:)=sum(zombiecountableTr((2*i-1):2*i,:),1);
end
pooledguides(pooledguides>1)=1;
pooledguides=[pooledguides;zombiecountableTr(71:end,:)];
writematrix(pooledguides,'\\helium\broad_clearylab\Users\Loic\thp1homedish2_zombie\properneighbors\pooledguidestableTER.csv')
pooledguides(36,:)=pooledguides(37,:)+pooledguides(38,:)+pooledguides(39,:);
pooledguides=pooledguides(1:36,:);
writematrix(pooledguides,'\\helium\broad_clearylab\Users\Loic\thp1homedish2_zombie\properneighbors\pooledguidestableTERjustonecontrol.csv')

