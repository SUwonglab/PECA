Group1=importdata('Sample1');
Group2=importdata('Sample2');
for grp=1:size(Group1,1)
	Node1_1{1,grp}=importdata(['../',Group1{grp,1},'/TFName.txt']);
	Node2_1{1,grp}=importdata(['../',Group1{grp,1},'/TGName.txt']);
	E=dlmread(['../',Group1{grp,1},'/TFTG_regulationScore.txt'],'\t');
	E(isnan(E))=0;
        fileID = fopen(['../../Input/',Group1{grp,1},'.txt']);
	C = textscan(fileID,'%s %f32');
	fclose(fileID);
	Symbol=C{1,1};
	G=C{1,2};
	[d1 f1]=ismember(Node1_1{1,grp},Symbol);
	[d2 f2]=ismember(Node2_1{1,grp},Symbol);
	Exp1=G(f1);
	Exp2=G(f2);	
	E=repmat(1./sqrt(Exp1+0.01),1,size(E,2)).*E.*repmat(1./sqrt(Exp2'+0.01),size(E,1),1);
	a=zscore(E);
	b=zscore(E')';
	a(a<0)=0;
	b(b<0)=0;
	E=a.*b;
	E1{1,grp}=E;
	E1_median(1,grp)=median(E(E>0));
	if exist('Group1_TF')~=1
		Group1_TF =Node1_1{1,grp};
		Group1_TG =Node2_1{1,grp};
	end
	Group1_TF=intersect(Group1_TF,Node1_1{1,grp});
	Group1_TG=intersect(Group1_TG,Node2_1{1,grp});
	Exp1_1{1,grp}=Exp1;
	Exp2_1{1,grp}=Exp2;
end
for grp=1:size(Group2,1)
	Node1_2{1,grp}=importdata(['../',Group2{grp,1},'/TFName.txt']);
	Node2_2{1,grp}=importdata(['../',Group2{grp,1},'/TGName.txt']);
	E=dlmread(['../',Group2{grp,1},'/TFTG_regulationScore.txt'],'\t');
	E(isnan(E))=0;
	fileID = fopen(['../../Input/',Group2{grp,1},'.txt']);
	C = textscan(fileID,'%s %f32');
	fclose(fileID);
	Symbol=C{1,1};
	G=C{1,2};
	[d1 f1]=ismember(Node1_2{1,grp},Symbol);
	[d2 f2]=ismember(Node2_2{1,grp},Symbol);
	Exp1=G(f1);
	Exp2=G(f2);	
	E=repmat(1./sqrt(Exp1+0.01),1,size(E,2)).*E.*repmat(1./sqrt(Exp2'+0.01),size(E,1),1);	
	a=zscore(E);
	b=zscore(E')';
	a(a<0)=0;
	b(b<0)=0;
	E=a.*b;
	E2{1,grp}=E;
	E2_median(1,grp)=median(E(E>0));
	if exist('Group2_TF')~=1
		Group2_TF =Node1_2{1,grp};
		Group2_TG =Node2_2{1,grp};
	end
	Group2_TF=intersect(Group2_TF,Node1_2{1,grp});
	Group2_TG=intersect(Group2_TG,Node2_2{1,grp});
	Exp1_2{1,grp}=Exp1;
	Exp2_2{1,grp}=Exp2;
end
E_mm=mean([E1_median E2_median]);
Node1=intersect(Group1_TF,Group2_TF);
Node2=intersect(Group1_TG,Group2_TG);
for grp=1:size(Group1,1)
	%E1{1,grp}=E1{1,grp}*E_mm/E1_median(1,grp);
	[d f1]=ismember(Node1,Node1_1{1,grp});
	[d f2]=ismember(Node2,Node2_1{1,grp});
	E1{1,grp}=E1{1,grp}(f1,f2);
	NetPool(grp,:,:)=E1{1,grp};
	ExpPool1(:,grp)=Exp1_1{1,grp}(f1);
	ExpPool2(:,grp)=Exp2_1{1,grp}(f2);	
end
for grp=1:size(Group2,1)
	%E2{1,grp}=E2{1,grp}*E_mm/E2_median(1,grp);
	[d f1]=ismember(Node1,Node1_2{1,grp});
	[d f2]=ismember(Node2,Node2_2{1,grp});
	E2{1,grp}=E2{1,grp}(f1,f2);
	NetPool(grp+size(Group1,1),:,:)=E2{1,grp};
	ExpPool1(:,grp+size(Group1,1))=Exp1_2{1,grp}(f1);
	ExpPool2(:,grp+size(Group1,1))=Exp2_2{1,grp}(f2);	
end
ExpPool1=double(ExpPool1);
ExpPool2=double(ExpPool2);
%%%%%%%%%%%%%%%%Edge
Net1_mean=reshape(mean(NetPool(1:size(Group1,1),:,:),1),length(Node1),length(Node2));
Net2_mean=reshape(mean(NetPool(1+size(Group1,1):end,:,:),1),length(Node1),length(Node2));
d1=max(ExpPool1')'>10;
d2=max(ExpPool2')'>10;
Node1=Node1(d1);
Node2=Node2(d2);
Net1_mean=Net1_mean(d1,d2);
Net2_mean=Net2_mean(d1,d2);
NetPool=NetPool(:,d1,d2);
ExpPool1=ExpPool1(d1,:);
ExpPool2=ExpPool2(d2,:);
load('../../Prior/TFTG_corr_organism.mat')
[d f2]=ismember(Node2,List);
[d f1]=ismember(Node1,TFName);
R2=R2(f1,f2);
Net1_mean=Net1_mean/(max(max(Net1_mean))+eps);
Net2_mean=Net2_mean/(max(max(Net2_mean))+eps);
id=find(abs(log2((Net1_mean+0.01)./(Net2_mean+0.01)))>log2(1.2));
NetPool_reshape=reshape(NetPool,size(Group1,1)+size(Group2,1),length(Node1)*length(Node2))';
p1=dif_test(NetPool_reshape(:,1:size(Group1,1)),NetPool_reshape(:,1+size(Group1,1):end));
p1(p1<10^(-16))=10^(-16);
p_reshape=-log10(reshape(p1,length(Node1),length(Node2)));
Net1=p_reshape.*(Net1_mean+0.01)./(Net2_mean+0.01).*Net1_mean/max(max(Net1_mean)).*(Net1_mean>Net2_mean);
Net2=p_reshape.*(Net2_mean+0.01)./(Net1_mean+0.01).*Net2_mean/max(max(Net2_mean)).*(Net2_mean>Net1_mean);
%%%%%%%%%%Node
q1=dif_test(ExpPool1(:,1:size(Group1,1)),ExpPool1(:,1+size(Group1,1):end));
q2=dif_test(ExpPool2(:,1:size(Group1,1)),ExpPool2(:,1+size(Group1,1):end));
q1(q1<10^(-16))=10^(-16);
q2(q2<10^(-16))=10^(-16);
q1=-log10(q1);
q2=-log10(q2);
%%%Net1 specific 
D1=(Net1>prctile(Net1(Net1>0),90)).*repmat(mean(ExpPool1(:,1:size(Group1,1)),2)>mean(ExpPool1(:,1+size(Group1,1):end),2),1,size(Net1,2));
[a b]=find(D1>0);
c=Net1(D1>0).*q1(a).*q2(b).*abs(R2(D1>0));
Net1_specific=[Node1(a) Node2(b) num2cell(R2(D1>0)) num2cell(q1(a)) num2cell(q2(b)) num2cell(p_reshape(D1>0)) num2cell((Net1_mean(D1>0)+0.01)./(Net2_mean(D1>0)+0.01)) num2cell(Net1_mean(D1>0)/max(max(Net1_mean)))];
d=R2(D1>0).*sign(mean(ExpPool2(b,1:size(Group1,1)),2)-mean(ExpPool2(b,1+size(Group1,1):end),2))>0;
Net1_specific=Net1_specific(d,:);
[d f]=sort(c(d),'descend');
Net1_specific=[Net1_specific(f,:) num2cell(d)];
Net1_specific=Net1_specific(abs(cell2mat(Net1_specific(:,3)))>0.1,:);
w1=double(q1').*(mean(ExpPool1(:,1:size(Group1,1)),2)>mean(ExpPool1(:,1+size(Group1,1):end),2))';
w2=double(q2');
[XNext, YNext,Obj]= FindSubNetwork_bi2(1,0.01,double(Net1),w1,w2);
d1=ismember(Net1_specific(:,1:2),[Node1(XNext>prctile(XNext,90));Node2(YNext>prctile(YNext,90))]);
Net1_specific_module=Net1_specific(d1(:,1).*d1(:,2)>0,:);
%%Net2 specific 
D1=(Net2>prctile(Net2(Net2>0),90)).*repmat(mean(ExpPool1(:,1:size(Group1,1)),2)<mean(ExpPool1(:,1+size(Group1,1):end),2),1,size(Net1,2));
[a b]=find(D1>0);
c=Net2(D1>0).*q1(a).*q2(b).*abs(R2(D1>0));
Net2_specific=[Node1(a) Node2(b) num2cell(R2(D1>0)) num2cell(q1(a)) num2cell(q2(b)) num2cell(p_reshape(D1>0)) num2cell((Net2_mean(D1>0)+1)./(Net1_mean(D1>0)+1)) num2cell(Net2_mean(D1>0)/max(max(Net2_mean)))];
d=R2(D1>0).*sign(mean(ExpPool2(b,1+size(Group1,1):end),2)-mean(ExpPool2(b,1:size(Group1,1)),2))>0;
Net2_specific=Net2_specific(d,:);
[d f]=sort(c(d),'descend');
Net2_specific=[Net2_specific(f,:) num2cell(d)];
Net2_specific=Net2_specific(abs(cell2mat(Net2_specific(:,3)))>0.1,:);
w1=double(q1').*(mean(ExpPool1(:,1:size(Group1,1)),2)'<mean(ExpPool1(:,1+size(Group1,1):end),2)');
w2=double(q2');
[XNext, YNext,Obj]= FindSubNetwork_bi2(1,0.01,double(Net2),w1,w2);
d1=ismember(Net2_specific(:,1:2),[Node1(XNext>prctile(XNext,90));Node2(YNext>prctile(YNext,90))]);
Net2_specific_module=Net2_specific(d1(:,1).*d1(:,2)>0,:);
%%%%%%%%%%%%%%%%common
Net3=sqrt(Net1_mean/max(max(Net1_mean)).*Net2_mean/max(max(Net2_mean)));
D1=Net3>max(prctile(Net3(Net3>0),95),0.01);
[a b]=find(D1>0);
c=Net1_mean(D1>0)/max(max(Net1_mean)).*Net2_mean(D1>0)/max(max(Net2_mean));
Net_common=[Node1(a) Node2(b) num2cell(R2(D1>0)) num2cell(Net1_mean(D1>0)/max(max(Net1_mean))) num2cell(Net2_mean(D1>0)/max(max(Net2_mean)))];
[d f]=sort(c,'descend');
Net_common=[Net_common(f,:) num2cell(d)];
Net_common=Net_common(abs(cell2mat(Net_common(:,3)))>0.1,:);
[XNext, YNext,Obj]= FindSubNetwork_bi1(2,double(Net3));
d1=ismember(Net_common(:,1:2),[Node1(XNext>prctile(XNext,95));Node2(YNext>prctile(YNext,95))]);
Net_common_module=Net_common(d1(:,1).*d1(:,2)>0,:);
%%%%%%%%%%%%%%%%write
Node_label=1*(mean(ExpPool2(:,1:size(Group1,1)),2)>mean(ExpPool2(:,1+size(Group1,1):end),2));
Node_label(Node_label==0)=2;
filename='Node_label.txt';
fid=fopen(filename,'wt');
for i=1:size(Node2,1)
	fprintf(fid, '%s\t',Node2{i,1});
	fprintf(fid, '%g\n',Node_label(i,1));
end
fclose(fid);
filename='Sample1_specific_network.txt';
fid=fopen(filename,'wt');
	fprintf(fid, '%s\t','TF');
	fprintf(fid, '%s\t','TG');
	fprintf(fid, '%s\t','corr');
	fprintf(fid, '%s\t','-log10P_TF');
	fprintf(fid, '%s\t','-log10P_TG');
	fprintf(fid, '%s\t','-log10P_Regulation');
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
	fprintf(fid, '%g\t',Net1_specific{i,6});
	fprintf(fid, '%g\t',Net1_specific{i,7});
	fprintf(fid, '%g\t',Net1_specific{i,8});
	fprintf(fid, '%g\n',Net1_specific{i,9});
%	fprintf(fid, '%s\n',Net1_specific{i,10});
end
fclose(fid);
filename='Sample2_specific_network.txt';
fid=fopen(filename,'wt');
	fprintf(fid, '%s\t','TF');
	fprintf(fid, '%s\t','TG');
	fprintf(fid, '%s\t','corr');
	fprintf(fid, '%s\t','-log10P_TF');
	fprintf(fid, '%s\t','-log10P_TG');
	fprintf(fid, '%s\t','-log10P_Regulation');
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
	fprintf(fid, '%g\t',Net2_specific{i,6});
	fprintf(fid, '%g\t',Net2_specific{i,7});
	fprintf(fid, '%g\t',Net2_specific{i,8});
	fprintf(fid, '%g\n',Net2_specific{i,9});
%	fprintf(fid, '%s\n',Net2_specific{i,10});
end
fclose(fid);
filename='Sample1_specific_module.txt';
fid=fopen(filename,'wt');
	fprintf(fid, '%s\t','TF');
	fprintf(fid, '%s\t','TG');
	fprintf(fid, '%s\t','corr');
	fprintf(fid, '%s\t','-log10P_TF');
	fprintf(fid, '%s\t','-log10P_TG');
	fprintf(fid, '%s\t','-log10P_Regulation');
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
	fprintf(fid, '%g\t',Net1_specific_module{i,6});
	fprintf(fid, '%g\t',Net1_specific_module{i,7});
	fprintf(fid, '%g\t',Net1_specific_module{i,8});
	fprintf(fid, '%g\n',Net1_specific_module{i,9});
%	fprintf(fid, '%s\n',Net1_specific_module{i,10});
end
fclose(fid);
filename='Sample2_specific_module.txt';
fid=fopen(filename,'wt');
	fprintf(fid, '%s\t','TF');
	fprintf(fid, '%s\t','TG');
	fprintf(fid, '%s\t','corr');
	fprintf(fid, '%s\t','-log10P_TF');
	fprintf(fid, '%s\t','-log10P_TG');
	fprintf(fid, '%s\t','-log10P_Regulation');
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
	fprintf(fid, '%g\t',Net2_specific_module{i,6});
	fprintf(fid, '%g\t',Net2_specific_module{i,7});
	fprintf(fid, '%g\t',Net2_specific_module{i,8});
	fprintf(fid, '%g\n',Net2_specific_module{i,9});
%	fprintf(fid, '%s\n',Net2_specific_module{i,10});
end
fclose(fid);
filename='Sample1_Sample2_common_network.txt';
fid=fopen(filename,'wt');
	fprintf(fid, '%s\t','TF');
	fprintf(fid, '%s\t','TG');
	fprintf(fid, '%s\t','corr');
	fprintf(fid, '%s\t','Activity_Sample1');
	fprintf(fid, '%s\t','Activity_Sample2');
	fprintf(fid, '%s\n','Activity');
for i=1:size(Net_common,1)
	fprintf(fid, '%s\t',Net_common{i,1});
	fprintf(fid, '%s\t',Net_common{i,2});
	fprintf(fid, '%g\t',Net_common{i,3});
	fprintf(fid, '%g\t',Net_common{i,4});
	fprintf(fid, '%g\t',Net_common{i,5});
	fprintf(fid, '%g\n',Net_common{i,6});
%	fprintf(fid, '%s\n',Net_common{i,7});
end
fclose(fid);
filename='Sample1_Sample2_common_module.txt';
fid=fopen(filename,'wt');
	fprintf(fid, '%s\t','TF');
	fprintf(fid, '%s\t','TG');
	fprintf(fid, '%s\t','corr');
	fprintf(fid, '%s\t','Activity_Sample1');
	fprintf(fid, '%s\t','Activity_Sample2');
	fprintf(fid, '%s\n','Activity');
for i=1:size(Net_common_module,1)
	fprintf(fid, '%s\t',Net_common_module{i,1});
	fprintf(fid, '%s\t',Net_common_module{i,2});
	fprintf(fid, '%g\t',Net_common_module{i,3});
	fprintf(fid, '%g\t',Net_common_module{i,4});
	fprintf(fid, '%g\t',Net_common_module{i,5});
	fprintf(fid, '%g\n',Net_common_module{i,6});
%	fprintf(fid, '%s\n',Net_common_module{i,4});
end
fclose(fid);

filename='Sample1_specific_TFnetwork.txt';
Net1_TFspecific=Net1_specific(ismember(Net1_specific(:,2),Node1),:);
fid=fopen(filename,'wt');
	fprintf(fid, '%s\t','TF');
	fprintf(fid, '%s\t','TG');
	fprintf(fid, '%s\t','corr');
	fprintf(fid, '%s\t','-log10P_TF');
	fprintf(fid, '%s\t','-log10P_TG');
	fprintf(fid, '%s\t','-log10P_Regulation');
	fprintf(fid, '%s\t','Fold');
	fprintf(fid, '%s\t','Activity');
	fprintf(fid, '%s\n','Score');
%	fprintf(fid, '%s\n','REs');
for i=1:size(Net1_TFspecific,1)
	fprintf(fid, '%s\t',Net1_TFspecific{i,1});
	fprintf(fid, '%s\t',Net1_TFspecific{i,2});
	fprintf(fid, '%g\t',Net1_TFspecific{i,3});
	fprintf(fid, '%g\t',Net1_TFspecific{i,4});
	fprintf(fid, '%g\t',Net1_TFspecific{i,5});
	fprintf(fid, '%g\t',Net1_TFspecific{i,6});
	fprintf(fid, '%g\t',Net1_TFspecific{i,7});
	fprintf(fid, '%g\t',Net1_TFspecific{i,8});
	fprintf(fid, '%g\n',Net1_TFspecific{i,9});
%	fprintf(fid, '%s\n',Net1_TFspecific{i,10});
end
fclose(fid);
filename='Sample2_specific_TFnetwork.txt';
Net2_TFspecific=Net2_specific(ismember(Net2_specific(:,2),Node1),:);
fid=fopen(filename,'wt');
	fprintf(fid, '%s\t','TF');
	fprintf(fid, '%s\t','TG');
	fprintf(fid, '%s\t','corr');
	fprintf(fid, '%s\t','-log10P_TF');
	fprintf(fid, '%s\t','-log10P_TG');
	fprintf(fid, '%s\t','-log10P_Regulation');
	fprintf(fid, '%s\t','Fold');
	fprintf(fid, '%s\t','Activity');
	fprintf(fid, '%s\n','Score');
%	fprintf(fid, '%s\n','REs');
for i=1:size(Net2_TFspecific,1)
	fprintf(fid, '%s\t',Net2_TFspecific{i,1});
	fprintf(fid, '%s\t',Net2_TFspecific{i,2});
	fprintf(fid, '%g\t',Net2_TFspecific{i,3});
	fprintf(fid, '%g\t',Net2_TFspecific{i,4});
	fprintf(fid, '%g\t',Net2_TFspecific{i,5});
	fprintf(fid, '%g\t',Net2_TFspecific{i,6});
	fprintf(fid, '%g\t',Net2_TFspecific{i,7});
	fprintf(fid, '%g\t',Net2_TFspecific{i,8});
	fprintf(fid, '%g\n',Net2_TFspecific{i,9});
%	fprintf(fid, '%s\n',Net2_TFspecific{i,10});
end
fclose(fid);
