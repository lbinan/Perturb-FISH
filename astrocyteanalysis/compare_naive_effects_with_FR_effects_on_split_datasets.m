split1=h5read('\\helium\broad_clearylab\Users\Loic\AsdRevstables\goodgenes\split1_count_data.h5ad','/X');
split1_perturbation_design=readtable('\\helium\broad_clearylab\Users\Loic\AsdRevstables\goodgenes\split1_perturbation_design.csv');
split1_perturbation_design=split1_perturbation_design(:,2:end);
split1_perturbation_design=table2array(split1_perturbation_design);
split1_perturbation_design=split1_perturbation_design';
split1=split1';
split1_perturbation_design(:,128)=((split1_perturbation_design(:,128)+split1_perturbation_design(:,129))>0)+0;
split1_perturbation_design=split1_perturbation_design(:,1:128);
for i=1:127
    effectssplit1(:,i)=mean(split1(split1_perturbation_design(:,i)==1,:),1)./mean(split1(split1_perturbation_design(:,128)==1,:),1);
    [H, psplit1(:,i)]=ttest2(split1(split1_perturbation_design(:,i)==1,:),split1(split1_perturbation_design(:,128)==1,:));
end
%%
split2=h5read('\\helium\broad_clearylab\Users\Loic\AsdRevstables\goodgenes\split2_count_data.h5ad','/X');
split2_perturbation_design=readtable('\\helium\broad_clearylab\Users\Loic\AsdRevstables\goodgenes\split2_perturbation_design.csv');
split2_perturbation_design=split2_perturbation_design(:,2:end);
split2_perturbation_design=table2array(split2_perturbation_design);
split2_perturbation_design=split2_perturbation_design';
split2=split2';
split2_perturbation_design(:,128)=((split2_perturbation_design(:,128)+split2_perturbation_design(:,129))>0)+0;
split2_perturbation_design=split2_perturbation_design(:,1:128);
for i=1:127
    effectssplit2(:,i)=mean(split2(split2_perturbation_design(:,i)==1,:),1)./mean(split2(split2_perturbation_design(:,128)==1,:),1);
    [H, psplit2(:,i)]=ttest2(split2(split2_perturbation_design(:,i)==1,:),split2(split2_perturbation_design(:,128)==1,:));
end
%%
LFC1=readtable('\\helium\broad_clearylab\Users\Loic\AsdRevstables\goodgenes\split1_LFCs.csv');
LFC1=LFC1(:,2:end);
LFC1=table2array(LFC1);
Q1=readtable('\\helium\broad_clearylab\Users\Loic\AsdRevstables\goodgenes\split1_qvals.csv');
Q1=Q1(:,2:end);
Q1=table2array(Q1);
LFC1=LFC1(:,1:127);
Q1=Q1(:,1:127);

LFC2=readtable('\\helium\broad_clearylab\Users\Loic\AsdRevstables\goodgenes\split2_LFCs.csv');
LFC2=LFC2(:,2:end);
LFC2=table2array(LFC2);
Q2=readtable('\\helium\broad_clearylab\Users\Loic\AsdRevstables\goodgenes\split2_qvals.csv');
Q2=Q2(:,2:end);
Q2=table2array(Q2);
LFC2=LFC2(:,1:127);
Q2=Q2(:,1:127);

LFC1_flat=reshape(LFC1,[127*277,1]);
LFC2_flat=reshape(LFC2,[127*277,1]);
Q1_flat=reshape(Q1,[127*277,1]);
Q2_flat=reshape(Q2,[127*277,1]);

effectssplit1_flat=log(reshape(effectssplit1,[127*277,1]));
effectssplit2_flat=log(reshape(effectssplit2,[127*277,1]));
psplit1_flat=reshape(psplit1,[127*277,1]);
psplit2_flat=reshape(psplit2,[127*277,1]);

subplot(2,2,1)
scatter(effectssplit1_flat,effectssplit2_flat,'.')
hold on 
scatter(effectssplit1_flat(psplit1_flat<0.1&psplit2_flat<0.1),effectssplit2_flat(psplit1_flat<0.1&psplit2_flat<0.1),'.')
hold off
xlabel('split1')
ylabel('split2')
C=corrcoef(effectssplit1_flat(psplit1_flat<0.1&psplit2_flat<0.1&~isinf(effectssplit1_flat)&~isinf(effectssplit2_flat)),effectssplit2_flat(psplit1_flat<0.1&psplit2_flat<0.1&~isinf(effectssplit1_flat)&~isinf(effectssplit2_flat)));
annotation('textbox',[0.15 0.9 0 0],'string',strcat('r=',num2str(C(1,2))),'FitBoxToText','on','EdgeColor','black')     

title('naive effects')

subplot(2,2,2)
scatter(effectssplit1_flat,LFC1_flat,'.')
hold on 
scatter(effectssplit1_flat(psplit1_flat<0.1&Q1_flat<0.1),LFC1_flat(psplit1_flat<0.1&Q1_flat<0.1),'.')
hold off
xlabel('naive')
ylabel('FR perturb')
C=corrcoef(effectssplit1_flat(psplit1_flat<0.1&Q1_flat<0.1&~isinf(effectssplit1_flat)),LFC1_flat(psplit1_flat<0.1&Q1_flat<0.1&~isinf(effectssplit1_flat)));
annotation('textbox',[0.6 0.9 0 0],'string',strcat('r=',num2str(C(1,2))),'FitBoxToText','on','EdgeColor','black')     
title('Split 1')

subplot(2,2,3)
scatter(effectssplit2_flat,LFC2_flat,'.')
hold on 
scatter(effectssplit2_flat(psplit2_flat<0.1&Q2_flat<0.1),LFC2_flat(psplit2_flat<0.1&Q2_flat<0.1),'.')
hold off
xlabel('naive')
ylabel('FR perturb')
C=corrcoef(effectssplit2_flat(psplit2_flat<0.1&Q2_flat<0.1),LFC2_flat(psplit2_flat<0.1&Q2_flat<0.1));
annotation('textbox',[0.15 0.45 0 0],'string',strcat('r=',num2str(C(1,2))),'FitBoxToText','on','EdgeColor','black')     

title('Split 2')

subplot(2,2,4)
scatter(LFC1_flat,LFC2_flat,'.')
hold on 
scatter(LFC1_flat(Q1_flat<0.1&Q2_flat<0.1),LFC2_flat(Q1_flat<0.1&Q2_flat<0.1),'.')
hold off
xlabel('split 1')
ylabel('split 2')
C=corrcoef(LFC1_flat(Q1_flat<0.1&Q2_flat<0.1),LFC2_flat(Q1_flat<0.1&Q2_flat<0.1));
annotation('textbox',[0.6 0.45 0 0],'string',strcat('r=',num2str(C(1,2))),'FitBoxToText','on','EdgeColor','black')     

title('FR perturb')
set(gcf,'color','w');

%% repeat on thp1s
split1=h5read('\\helium\broad_clearylab\Users\Loic\dataforbrian\zombie_sample1\split1_count_data.h5ad','/X');
split1_perturbation_design=readtable('\\helium\broad_clearylab\Users\Loic\dataforbrian\zombie_sample1\split1_perturbation_design.csv');
split1_perturbation_design=split1_perturbation_design(:,2:end);
split1_perturbation_design=table2array(split1_perturbation_design);
split1_perturbation_design=split1_perturbation_design';
split1=split1';
for i=1:34
    effectssplit1(:,i)=mean(split1(split1_perturbation_design(:,i)==1,:),1)./mean(split1(split1_perturbation_design(:,35)==1,:),1);
    [H, psplit1(:,i)]=ttest2(split1(split1_perturbation_design(:,i)==1,:),split1(split1_perturbation_design(:,35)==1,:));
end
%%
split2=h5read('\\helium\broad_clearylab\Users\Loic\dataforbrian\zombie_sample1\split2_count_data.h5ad','/X');
split2_perturbation_design=readtable('\\helium\broad_clearylab\Users\Loic\dataforbrian\zombie_sample1\split2_perturbation_design.csv');
split2_perturbation_design=split2_perturbation_design(:,2:end);
split2_perturbation_design=table2array(split2_perturbation_design);
split2_perturbation_design=split2_perturbation_design';
split2=split2';
for i=1:34
    effectssplit2(:,i)=mean(split2(split2_perturbation_design(:,i)==1,:),1)./mean(split2(split2_perturbation_design(:,35)==1,:),1);
    [H, psplit2(:,i)]=ttest2(split2(split2_perturbation_design(:,i)==1,:),split2(split2_perturbation_design(:,35)==1,:));
end
%%
LFC1=readtable('\\helium\broad_clearylab\Users\Loic\dataforbrian\zombie_sample1\split1_LFCs.csv');
LFC1=LFC1(:,2:end);
LFC1=table2array(LFC1);
Q1=readtable('\\helium\broad_clearylab\Users\Loic\dataforbrian\zombie_sample1\split1_qvals.csv');
Q1=Q1(:,2:end);
Q1=table2array(Q1);
LFC1=LFC1(:,1:34);
Q1=Q1(:,1:34);

LFC2=readtable('\\helium\broad_clearylab\Users\Loic\dataforbrian\zombie_sample1\split2_LFCs.csv');
LFC2=LFC2(:,2:end);
LFC2=table2array(LFC2);
Q2=readtable('\\helium\broad_clearylab\Users\Loic\dataforbrian\zombie_sample1\split2_qvals.csv');
Q2=Q2(:,2:end);
Q2=table2array(Q2);
LFC2=LFC2(:,1:34);
Q2=Q2(:,1:34);

LFC1_flat=reshape(LFC1,[128*34,1]);
LFC2_flat=reshape(LFC2,[128*34,1]);
Q1_flat=reshape(Q1,[128*34,1]);
Q2_flat=reshape(Q2,[128*34,1]);

effectssplit1_flat=log(reshape(effectssplit1,[128*34,1]));
effectssplit2_flat=log(reshape(effectssplit2,[128*34,1]));
psplit1_flat=reshape(psplit1,[128*34,1]);
psplit2_flat=reshape(psplit2,[128*34,1]);

subplot(2,2,1)
scatter(effectssplit1_flat,effectssplit2_flat,'.')
hold on 
scatter(effectssplit1_flat(psplit1_flat<0.1&psplit2_flat<0.1),effectssplit2_flat(psplit1_flat<0.1&psplit2_flat<0.1),'.')
hold off
xlabel('split1')
ylabel('split2')
C=corrcoef(effectssplit1_flat(psplit1_flat<0.1&psplit2_flat<0.1&~isinf(effectssplit1_flat)&~isinf(effectssplit2_flat)),effectssplit2_flat(psplit1_flat<0.1&psplit2_flat<0.1&~isinf(effectssplit1_flat)&~isinf(effectssplit2_flat)));
annotation('textbox',[0.15 0.9 0 0],'string',strcat('r=',num2str(C(1,2))),'FitBoxToText','on','EdgeColor','black')     

title('naive effects')

subplot(2,2,2)
scatter(effectssplit1_flat,LFC1_flat,'.')
hold on 
scatter(effectssplit1_flat(psplit1_flat<0.1&Q1_flat<0.1),LFC1_flat(psplit1_flat<0.1&Q1_flat<0.1),'.')
hold off
xlabel('naive')
ylabel('FR perturb')
C=corrcoef(effectssplit1_flat(psplit1_flat<0.1&Q1_flat<0.1&~isinf(effectssplit1_flat)),LFC1_flat(psplit1_flat<0.1&Q1_flat<0.1&~isinf(effectssplit1_flat)));
annotation('textbox',[0.6 0.9 0 0],'string',strcat('r=',num2str(C(1,2))),'FitBoxToText','on','EdgeColor','black')     
title('Split 1')

subplot(2,2,3)
scatter(effectssplit2_flat,LFC2_flat,'.')
hold on 
scatter(effectssplit2_flat(psplit2_flat<0.1&Q2_flat<0.1),LFC2_flat(psplit2_flat<0.1&Q2_flat<0.1),'.')
hold off
xlabel('naive')
ylabel('FR perturb')
C=corrcoef(effectssplit2_flat(psplit2_flat<0.1&Q2_flat<0.1),LFC2_flat(psplit2_flat<0.1&Q2_flat<0.1));
annotation('textbox',[0.15 0.45 0 0],'string',strcat('r=',num2str(C(1,2))),'FitBoxToText','on','EdgeColor','black')     

title('Split 2')

subplot(2,2,4)
scatter(LFC1_flat,LFC2_flat,'.')
hold on 
scatter(LFC1_flat(Q1_flat<0.1&Q2_flat<0.1),LFC2_flat(Q1_flat<0.1&Q2_flat<0.1),'.')
hold off
xlabel('split 1')
ylabel('split 2')
C=corrcoef(LFC1_flat(Q1_flat<0.1&Q2_flat<0.1),LFC2_flat(Q1_flat<0.1&Q2_flat<0.1));
annotation('textbox',[0.6 0.45 0 0],'string',strcat('r=',num2str(C(1,2))),'FitBoxToText','on','EdgeColor','black')     

title('FR perturb')
set(gcf,'color','w');

