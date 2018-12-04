Sample{1,1}='Sample1';
Sample{2,1}='Sample2';
%%%Sample1
Node1_1=importdata(['../',Sample{1,1},'/TFName.txt']);
Node2_1=importdata(['../',Sample{1,1},'/TGName.txt']);
E=dlmread(['../',Sample{1,1},'/TFTG_regulationScore.txt'],'\t');
E(isnan(E))=0;
a=zscore(E);
b=zscore(E')';
a(a<0)=0;
b(b<0)=0;
E=a.*b;
E1{1,1}=E;
E_median(1,1)=median(E(E>0));
fileID = fopen(['../../Input/',Sample{1,1},'.txt']);
C = textscan(fileID,'%s %f32');
fclose(fileID);
Symbol=C{1,1};
G=C{1,2};
[d1 f1]=ismember(Node1_1,Symbol);
[d2 f2]=ismember(Node2_1,Symbol);
Exp1_1=G(f1);
Exp2_1=G(f2);
%%%Sample2
Node1_2=importdata(['../',Sample{2,1},'/TFName.txt']);
Node2_2=importdata(['../',Sample{2,1},'/TGName.txt']);	
E=dlmread(['../',Sample{2,1},'/TFTG_regulationScore.txt'],'\t');
E(isnan(E))=0;
a=zscore(E);
b=zscore(E')';
a(a<0)=0;
b(b<0)=0;
E=a.*b;
E1{1,2}=E;
E_median(2,1)=median(E(E>0));
fileID = fopen(['../../Input/',Sample{2,1},'.txt']);
C = textscan(fileID,'%s %f32');
fclose(fileID);
Symbol=C{1,1};
G=C{1,2};
[d1 f1]=ismember(Node1_2,Symbol);
[d2 f2]=ismember(Node2_2,Symbol);
Exp1_2=G(f1);
Exp2_2=G(f2);
%%%%%
E_mm=mean(E_median);
for i=1:2
E1{1,i}=E1{1,i}*E_mm/E_median(i);
end
Node1=intersect(Node1_1,Node1_2);
Node2=intersect(Node2_1,Node2_2);
[d f1]=ismember(Node1,Node1_1);
[d f2]=ismember(Node2,Node2_1);
[d f3]=ismember(Node1,Node1_2);
[d f4]=ismember(Node2,Node2_2);
E1{1,1}=E1{1,1}(f1,f2);
E1{1,2}=E1{1,2}(f3,f4);
ExpPool1=[Exp1_1(f1) Exp1_2(f3)];
ExpPool2=[Exp2_1(f2) Exp2_2(f4)];
%%%%%%%%%
d1=max(ExpPool1')'>10;
d2=max(ExpPool2')'>10;
Node1=Node1(d1);
Node2=Node2(d2);
E1{1,1}=E1{1,1}(d1,d2);
E1{1,2}=E1{1,2}(d1,d2);
ExpPool1=ExpPool1(d1,:);
ExpPool2=ExpPool2(d2,:);
load('../../Prior/TFTG_corr_organism.mat')
[d f2]=ismember(Node2,List);
[d f1]=ismember(Node1,TFName);
R2=R2(f1,f2);
Net1=(E1{1,1}>E1{1,2}).*E1{1,1}.*(E1{1,1}+1)./(E1{1,2}+1);
Net2=(E1{1,2}>E1{1,1}).*E1{1,2}.*(E1{1,2}+1)./(E1{1,1}+1);
Net1=Net1/max(max(Net1));
Net2=Net2/max(max(Net2));
fileID = fopen(['../',Sample{1,1},'/',Sample{1,1},'_network.txt']);
C = textscan(fileID,'%s %s %f32 %s','HeaderLines', 1);
fclose(fileID);
Net1_all=[C{1,1} C{1,2} num2cell(C{1,3}) C{1,4}];
fileID = fopen(['../',Sample{2,1},'/',Sample{2,1},'_network.txt']);
C = textscan(fileID,'%s %s %f32 %s','HeaderLines', 1);
fclose(fileID);
Net2_all=[C{1,1} C{1,2} num2cell(C{1,3}) C{1,4}];
%%%Net1 specific 
D1=(Net1>prctile(Net1(Net1>0),90)).*repmat(ExpPool1(:,1)>ExpPool1(:,2),1,size(Net1,2));
[a b]=find(D1>0);
c=Net1(D1>0).*abs(R2(D1>0));
Net1_specific=[Node1(a) Node2(b) num2cell(R2(D1>0)) num2cell((Net1(D1>0)+0.01)./(Net2(D1>0)+0.01)) num2cell(Net1(D1>0)/max(max(Net1)))];
d=R2(D1>0).*sign(ExpPool2(b,1)-ExpPool2(b,2))>0;
Net1_specific=Net1_specific(d,:);
[d f]=sort(c(d),'descend');
Net1_specific=[Net1_specific(f,:) num2cell(d)];
Net1_specific=Net1_specific(abs(cell2mat(Net1_specific(:,3)))>0.1,:);
[XNext, YNext,Obj]= FindSubNetwork_bi1(1,Net1);
d1=ismember(Net1_specific(:,1:2),[Node1(XNext>prctile(XNext,90));Node2(YNext>prctile(YNext,90))]);
Net1_specific_module=Net1_specific(d1(:,1).*d1(:,2)>0,:);
%%Net2 specific 
D1=(Net2>prctile(Net1(Net1>0),90)).*repmat(ExpPool1(:,2)>ExpPool1(:,1),1,size(Net2,2));
[a b]=find(D1>0);
c=Net2(D1>0).*abs(R2(D1>0));
Net2_specific=[Node1(a) Node2(b) num2cell(R2(D1>0)) num2cell((Net2(D1>0)+0.01)./(Net1(D1>0)+0.01)) num2cell(Net2(D1>0)/max(max(Net2)))];
d=R2(D1>0).*sign(ExpPool2(b,2)-ExpPool2(b,1))>0;
Net2_specific=Net2_specific(d,:);
[d f]=sort(c(d),'descend');
Net2_specific=[Net2_specific(f,:) num2cell(d)];
Net2_specific=Net2_specific(abs(cell2mat(Net2_specific(:,3)))>0.1,:);
[XNext, YNext,Obj]= FindSubNetwork_bi1(1,Net2);
d1=ismember(Net2_specific(:,1:2),[Node1(XNext>prctile(XNext,90));Node2(YNext>prctile(YNext,90))]);
Net2_specific_module=Net2_specific(d1(:,1).*d1(:,2)>0,:);
%%%%%%%%%%%%%%%%common
Net3=E1{1,1}+E1{1,2};
D1=Net3>prctile(Net3(Net3>0),95);
[a b]=find(D1>0);
c=E1{1,1}(D1>0).*E1{1,2}(D1>0)/(max(max(E1{1,1}))*max(max(E1{1,2})));
Net_common=[Node1(a) Node2(b) num2cell(R2(D1>0)) num2cell(E1{1,1}(D1>0)/max(max(E1{1,1}))) num2cell(E1{1,2}(D1>0)/max(max(E1{1,2})))];
[d f]=sort(c,'descend');
Net_common=[Net_common(f,:) num2cell(d)];
Net_common=Net_common(abs(cell2mat(Net_common(:,3)))>0.1,:);
[XNext, YNext,Obj]= FindSubNetwork_bi1(1,Net3);
d1=ismember(Net_common(:,1:2),[Node1(XNext>prctile(XNext,95));Node2(YNext>prctile(YNext,95))]);
Net_common_module=Net_common(d1(:,1).*d1(:,2)>0,:);
%%%%%%%%%%%%%write
Node_label=1*(ExpPool2(:,1)>ExpPool2(:,2));
Node_label(Node_label==0)=2;
filename='Node_label.txt';
fid=fopen(filename,'wt');
for i=1:size(Node2,1)
	fprintf(fid, '%s\t',Node2{i,1});
	fprintf(fid, '%g\n',Node_label(i,1));
end
fclose(fid);

filename=[Sample{1,1},'_specific_network.txt'];
fid=fopen(filename,'wt');
	fprintf(fid, '%s\t','TF');
	fprintf(fid, '%s\t','TG');
	fprintf(fid, '%s\t','corr');
	fprintf(fid, '%s\t','Fold');
	fprintf(fid, '%s\t','Activity');
	fprintf(fid, '%s\n','Score');
for i=1:size(Net1_specific,1)
	fprintf(fid, '%s\t',Net1_specific{i,1});
	fprintf(fid, '%s\t',Net1_specific{i,2});
	fprintf(fid, '%g\t',Net1_specific{i,3});
	fprintf(fid, '%g\t',Net1_specific{i,4});
	fprintf(fid, '%g\t',Net1_specific{i,5});
	fprintf(fid, '%g\n',Net1_specific{i,6});
end
fclose(fid);
filename=[Sample{2,1},'_specific_network.txt'];
fid=fopen(filename,'wt');
	fprintf(fid, '%s\t','TF');
	fprintf(fid, '%s\t','TG');
	fprintf(fid, '%s\t','corr');
	fprintf(fid, '%s\t','Fold');
	fprintf(fid, '%s\t','Activity');
	fprintf(fid, '%s\n','Score');
for i=1:size(Net2_specific,1)
	fprintf(fid, '%s\t',Net2_specific{i,1});
	fprintf(fid, '%s\t',Net2_specific{i,2});
	fprintf(fid, '%g\t',Net2_specific{i,3});
	fprintf(fid, '%g\t',Net2_specific{i,4});
	fprintf(fid, '%g\t',Net2_specific{i,5});
	fprintf(fid, '%g\n',Net2_specific{i,6});
end
fclose(fid);
filename=[Sample{1,1},'_specific_module.txt'];
fid=fopen(filename,'wt');
	fprintf(fid, '%s\t','TF');
	fprintf(fid, '%s\t','TG');
	fprintf(fid, '%s\t','corr');
	fprintf(fid, '%s\t','Fold');
	fprintf(fid, '%s\t','Activity');
	fprintf(fid, '%s\n','Score');
for i=1:size(Net1_specific_module,1)
	fprintf(fid, '%s\t',Net1_specific_module{i,1});
	fprintf(fid, '%s\t',Net1_specific_module{i,2});
	fprintf(fid, '%g\t',Net1_specific_module{i,3});
	fprintf(fid, '%g\t',Net1_specific_module{i,4});
	fprintf(fid, '%g\t',Net1_specific_module{i,5});
	fprintf(fid, '%g\n',Net1_specific_module{i,6});
end
fclose(fid);
filename=[Sample{2,1},'_specific_module.txt'];
fid=fopen(filename,'wt');
	fprintf(fid, '%s\t','TF');
	fprintf(fid, '%s\t','TG');
	fprintf(fid, '%s\t','corr');
	fprintf(fid, '%s\t','Fold');
	fprintf(fid, '%s\t','Activity');
	fprintf(fid, '%s\n','Score');
for i=1:size(Net2_specific_module,1)
	fprintf(fid, '%s\t',Net2_specific_module{i,1});
	fprintf(fid, '%s\t',Net2_specific_module{i,2});
	fprintf(fid, '%g\t',Net2_specific_module{i,3});
	fprintf(fid, '%g\t',Net2_specific_module{i,4});
	fprintf(fid, '%g\t',Net2_specific_module{i,5});
	fprintf(fid, '%g\n',Net2_specific_module{i,6});
end
fclose(fid);
filename=[Sample{1,1},'_',Sample{2,1},'_common_network.txt'];
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
end
fclose(fid);
filename=[Sample{1,1},'_',Sample{2,1},'_common_module.txt'];
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
end
fclose(fid);
