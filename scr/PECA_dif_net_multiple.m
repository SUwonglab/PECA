Group1=importdata('Sample1');
Group2=importdata('Sample2');
for grp=1:size(Group1,1)
	Node1_1{1,grp}=importdata(['../',Group1{grp,1},'/TFName.txt']);
	Node2_1{1,grp}=importdata(['../',Group1{grp,1},'/TGName.txt']);
	E=dlmread(['../',Group1{grp,1},'/TFTG_regulationScore.txt'],'\t');
	E(isnan(E))=0;
	E1{1,grp}=E;
	E1_median(1,grp)=median(E(E>0));
	if exist('Group1_TF')~=1
		Group1_TF =Node1_1{1,grp};
		Group1_TG =Node2_1{1,grp};
	end
	Group1_TF=intersect(Group1_TF,Node1_1{1,grp});
	Group1_TG=intersect(Group1_TG,Node2_1{1,grp});
end
for grp=1:size(Group2,1)
	Node1_2{1,grp}=importdata(['../',Group2{grp,1},'/TFName.txt']);
	Node2_2{1,grp}=importdata(['../',Group2{grp,1},'/TGName.txt']);
	E=dlmread(['../',Group2{grp,1},'/TFTG_regulationScore.txt'],'\t');
	E(isnan(E))=0;
	E2{1,grp}=E;
	E2_median(1,grp)=median(E(E>0));
	if exist('Group2_TF')~=1
		Group2_TF =Node1_2{1,grp};
		Group2_TG =Node2_2{1,grp};
	end
	Group2_TF=intersect(Group2_TF,Node1_2{1,grp});
	Group2_TG=intersect(Group2_TG,Node2_2{1,grp});
end
E_mm=mean([E1_median E2_median]);
Node1=intersect(Group1_TF,Group2_TF);
Node2=intersect(Group1_TG,Group2_TG);
for grp=1:size(Group1,1)
	E1{1,grp}=E1{1,grp}*E_mm/E1_median(1,grp);
	[d f1]=ismember(Node1_1{1,grp},Node1);
	[d f2]=ismember(Node2_1{1,grp},Node2);
	E1{1,grp}=E1{1,grp}(f1,f2);
	NetPool(grp,:,:)=E1{1,grp};
end
for grp=1:size(Group2,1)
	E2{1,grp}=E2{1,grp}*E_mm/E2_median(1,grp);
	[d f1]=ismember(Node1_2{1,grp},Node1);
	[d f2]=ismember(Node2_2{1,grp},Node2);
	E2{1,grp}=E2{1,grp}(f1,f2);
	NetPool(grp+size(Group1,1),:,:)=E2{1,grp};
end
%%%%%%%%%%%%%%%%
Net1_mean=reshape(mean(NetPool(1:size(Group1,1),:,:),1),length(Node1),length(Node2));
Net2_mean=reshape(mean(NetPool(1+size(Group1,1):end,:,:),1),length(Node1),length(Node2));
id=find(abs(log2((Net1_mean+1)./(Net2_mean+1)))>log2(1.2));
NetPool_reshape=reshape(NetPool,size(Group1,1)+size(Group2,1),length(Node1)*length(Node2))';
p1=dif_test(NetPool_reshape(:,1:size(Group1,1)),NetPool_reshape(:,1+size(Group1,1):end));
p1(p1<10^(-16))=10^(-16);
p_reshape=-log10(reshape(p1,length(Node1),length(Node2)));
Net1=p_reshape.*(Net1_mean+1)./(Net2_mean+1).*Net1_mean/max(max(Net1_mean)).*(Net1_mean>Net2_mean).*(p_reshape>-log10(0.05));
Net2=p_reshape.*(Net2_mean+1)./(Net1_mean+1).*Net2_mean/max(max(Net2_mean)).*(Net2_mean>Net1_mean).*(p_reshape>-log10(0.05));
%%%Net1 specific 
D1=Net1>prctile(Net1(Net1>0),95);
[a b]=find(D1>0);
c=Net1(D1>0);
Net1_specific=[Node1(a) Node2(b) num2cell(p_reshape(D1>0)) num2cell((Net1_mean(D1>0)+1)./(Net2_mean(D1>0)+1)) num2cell(Net1_mean(D1>0)/max(max(Net1_mean)))];
[d f]=sort(c,'descend');
Net1_specific=[Net1_specific(f,:) num2cell(d)];
[XNext, YNext,Obj]= FindSubNetwork_bi1(1,Net1);
d1=ismember(Net1_specific(:,1),Node1(XNext>prctile(XNext,95)));
d1(:,2)=ismember(Net1_specific(:,2),[Node1(XNext>prctile(XNext,95));Node2(YNext>prctile(YNext,95))]);
Net1_specific_module=Net1_specific(d1(:,1).*d1(:,2)>0,:);
%%Net2 specific 
D1=Net2>prctile(Net2(Net2>0),95);
[a b]=find(D1>0);
c=Net2(D1>0);
Net2_specific=[Node1(a) Node2(b) num2cell(p_reshape(D1>0)) num2cell((Net2_mean(D1>0)+1)./(Net1_mean(D1>0)+1)) num2cell(Net2_mean(D1>0)/max(max(Net2_mean)))];
[d f]=sort(c,'descend');
Net2_specific=[Net2_specific(f,:) num2cell(d)];
[XNext, YNext,Obj]= FindSubNetwork_bi1(1,Net2);
d1=ismember(Net2_specific(:,1),Node1(XNext>prctile(XNext,95)));
d1(:,2)=ismember(Net2_specific(:,2),[Node1(XNext>prctile(XNext,95));Node2(YNext>prctile(YNext,95))]);
Net2_specific_module=Net2_specific(d1(:,1).*d1(:,2)>0,:);
%%%%write
filename='Sample1_specific_network.txt';
fid=fopen(filename,'wt');
	fprintf(fid, '%s\t','TF');
	fprintf(fid, '%s\t','TG');
	fprintf(fid, '%s\t','-log10P');
	fprintf(fid, '%s\t','Fold');
	fprintf(fid, '%s\t','Activity');
	fprintf(fid, '%s\n','Score');
%	fprintf(fid, '%s\n','REs');
for i=1:size(Net1_specific,1)
	fprintf(fid, '%s\t',Net1_specific{i,1});
	fprintf(fid, '%s\t',Net1_specific{i,2});
	fprintf(fid, '%g\t',Net1_specific{i,3});
	fprintf(fid, '%g\t',Net1_specific{i,4});
	fprintf(fid, '%g\t',Net1_specific{i,5});
	fprintf(fid, '%g\n',Net1_specific{i,6});
%	fprintf(fid, '%s\n',Net1_specific{i,7});
end
fclose(fid);
filename='Sample2_specific_network.txt';
fid=fopen(filename,'wt');
	fprintf(fid, '%s\t','TF');
	fprintf(fid, '%s\t','TG');
	fprintf(fid, '%s\t','-log10P');
	fprintf(fid, '%s\t','Fold');
	fprintf(fid, '%s\t','Activity');
	fprintf(fid, '%s\n','Score');
%	fprintf(fid, '%s\n','REs');
for i=1:size(Net2_specific,1)
	fprintf(fid, '%s\t',Net2_specific{i,1});
	fprintf(fid, '%s\t',Net2_specific{i,2});
	fprintf(fid, '%g\t',Net2_specific{i,3});
	fprintf(fid, '%g\t',Net2_specific{i,4});
	fprintf(fid, '%g\t',Net2_specific{i,5});
	fprintf(fid, '%g\n',Net2_specific{i,6});
%	fprintf(fid, '%s\n',Net2_specific{i,7});
end
fclose(fid);
filename='Sample1_specific_module.txt';
fid=fopen(filename,'wt');
	fprintf(fid, '%s\t','TF');
	fprintf(fid, '%s\t','TG');
	fprintf(fid, '%s\t','-log10P');
	fprintf(fid, '%s\t','Fold');
	fprintf(fid, '%s\t','Activity');
	fprintf(fid, '%s\n','Score');
%	fprintf(fid, '%s\n','REs');
for i=1:size(Net1_specific_module,1)
	fprintf(fid, '%s\t',Net1_specific_module{i,1});
	fprintf(fid, '%s\t',Net1_specific_module{i,2});
	fprintf(fid, '%g\t',Net1_specific_module{i,3});
	fprintf(fid, '%g\t',Net1_specific_module{i,4});
	fprintf(fid, '%g\t',Net1_specific_module{i,5});
	fprintf(fid, '%g\n',Net1_specific_module{i,6});
%	fprintf(fid, '%s\n',Net1_specific_module{i,7});
end
fclose(fid);
filename='Sample2_specific_module.txt';
fid=fopen(filename,'wt');
	fprintf(fid, '%s\t','TF');
	fprintf(fid, '%s\t','TG');
	fprintf(fid, '%s\t','-log10P');
	fprintf(fid, '%s\t','Fold');
	fprintf(fid, '%s\t','Activity');
	fprintf(fid, '%s\n','Score');
%	fprintf(fid, '%s\n','REs');
for i=1:size(Net2_specific_module,1)
	fprintf(fid, '%s\t',Net2_specific_module{i,1});
	fprintf(fid, '%s\t',Net2_specific_module{i,2});
	fprintf(fid, '%g\t',Net2_specific_module{i,3});
	fprintf(fid, '%g\t',Net2_specific_module{i,4});
	fprintf(fid, '%g\t',Net2_specific_module{i,5});
	fprintf(fid, '%g\n',Net2_specific_module{i,6});
%	fprintf(fid, '%s\n',Net2_specific_module{i,7});
end
fclose(fid);
