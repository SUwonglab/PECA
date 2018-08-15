Node1_1=importdata('../Sample1/TFName.txt');
Node2_1=importdata('../Sample1/TGName.txt');
E=dlmread(['../Sample1/TFTG_regulationScore.txt'],'\t');
E(isnan(E))=0;
E1{1,1}=E;
E_median(1,1)=median(E(E>0));
Node1_2=importdata('../Sample2/TFName.txt');
Node2_2=importdata('../Sample2/TGName.txt');
E=dlmread(['../Sample2/TFTG_regulationScore.txt'],'\t');
E(isnan(E))=0;
E1{1,2}=E;
E_median(2,1)=median(E(E>0));
E_mm=mean(E_median);
for i=1:2
E1{1,i}=E1{1,i}*E_mm/E_median(i);
end
Node1=intersect(Node1_1,Node1_2);
Node2=intersect(Node2_1,Node2_2);
E1{1,1}=E1{1,1}(ismember(Node1_1,Node1),ismember(Node2_1,Node2));
E1{1,2}=E1{1,2}(ismember(Node1_2,Node1),ismember(Node2_2,Node2));
Net1=sqrt(E1{1,1}.*(E1{1,1}+1)./(E1{1,2}+1));
Net2=sqrt(E1{1,2}.*(E1{1,2}+1)./(E1{1,1}+1));
fileID = fopen('../Sample1/Sample1_network.txt');
C = textscan(fileID,'%s %s %f32 %s','HeaderLines', 1);
fclose(fileID);
Net1_all=[C{1,1} C{1,2} num2cell(C{1,3}) C{1,4}];
fileID = fopen('../Sample2/Sample2_network.txt');
C = textscan(fileID,'%s %s %f32 %s','HeaderLines', 1);
fclose(fileID);
Net2_all=[C{1,1} C{1,2} num2cell(C{1,3}) C{1,4}];
%%%Net1 specific 
D1=Net1>prctile(Net1(Net1>0),99);
[a b]=find(D1>0);
c=Net1(D1>0);
Net1_specific=[Node1(a) Node2(b)];
[d f]=sort(c,'descend');
Net1_specific=[Net1_specific(f,:) num2cell(d)];
[XNext, YNext,Obj]= FindSubNetwork_bi1(1,Net1);
d1=ismember(Net1_specific(:,1),Node1(XNext>prctile(XNext,99)));
d1(:,2)=ismember(Net1_specific(:,2),[Node1(XNext>prctile(XNext,99));Node2(YNext>prctile(YNext,99))]);
Net1_specific_module=Net1_specific(d1(:,1).*d1(:,2)>0,:);
[d f]=ismember(Net1_all(:,1:2),Node2);
[d f1]=ismember(Net1_specific(:,1:2),Node2);
[d f2]=ismember(Net1_specific_module(:,1:2),Node2);
[dd ff]=ismember(f1,f,'rows');
Net1_specific=[Net1_specific(dd==1,:) Net1_all(ff(dd==1),4)];
[dd ff]=ismember(f2,f,'rows');
Net1_specific_module=[Net1_specific_module(dd==1,:) Net1_all(ff(dd==1),4)];
%%Net2 specific 
D1=Net2>prctile(Net2(Net2>0),99);
[a b]=find(D1>0);
c=Net2(D1>0);
Net2_specific=[Node1(a) Node2(b)];
[d f]=sort(c,'descend');
Net2_specific=[Net2_specific(f,:) num2cell(d)];
[XNext, YNext,Obj]= FindSubNetwork_bi1(1,Net2);
d1=ismember(Net2_specific(:,1),Node1(XNext>prctile(XNext,99)));
d1(:,2)=ismember(Net2_specific(:,2),[Node1(XNext>prctile(XNext,99));Node2(YNext>prctile(YNext,99))]);
Net2_specific_module=Net2_specific(d1(:,1).*d1(:,2)>0,:);
[d f]=ismember(Net2_all(:,1:2),Node2);
[d f1]=ismember(Net2_specific(:,1:2),Node2);
[d f2]=ismember(Net2_specific_module(:,1:2),Node2);
[dd ff]=ismember(f1,f,'rows');
Net2_specific=[Net2_specific(dd==1,:) Net2_all(ff(dd==1),4)];
[dd ff]=ismember(f2,f,'rows');
Net2_specific_module=[Net2_specific_module(dd==1,:) Net2_all(ff(dd==1),4)];
%%%%write
filename='Sample1_specific_network.txt';
fid=fopen(filename,'wt');
	fprintf(fid, '%s\t','TF');
	fprintf(fid, '%s\t','TG');
	fprintf(fid, '%s\t','Score');
	fprintf(fid, '%s\n','REs');
for i=1:size(Net1_specific,1)
	fprintf(fid, '%s\t',Net1_specific{i,1});
	fprintf(fid, '%s\t',Net1_specific{i,2});
	fprintf(fid, '%g\t',Net1_specific{i,3});
	fprintf(fid, '%s\n',Net1_specific{i,4});
end
fclose(fid);
filename='Sample2_specific_network.txt';
fid=fopen(filename,'wt');
	fprintf(fid, '%s\t','TF');
	fprintf(fid, '%s\t','TG');
	fprintf(fid, '%s\t','Score');
	fprintf(fid, '%s\n','REs');
for i=1:size(Net2_specific,1)
	fprintf(fid, '%s\t',Net2_specific{i,1});
	fprintf(fid, '%s\t',Net2_specific{i,2});
	fprintf(fid, '%g\t',Net2_specific{i,3});
	fprintf(fid, '%s\n',Net2_specific{i,4});
end
fclose(fid);
filename='Sample1_specific_module.txt';
fid=fopen(filename,'wt');
	fprintf(fid, '%s\t','TF');
	fprintf(fid, '%s\t','TG');
	fprintf(fid, '%s\t','Score');
	fprintf(fid, '%s\n','REs');
for i=1:size(Net1_specific_module,1)
	fprintf(fid, '%s\t',Net1_specific_module{i,1});
	fprintf(fid, '%s\t',Net1_specific_module{i,2});
	fprintf(fid, '%g\t',Net1_specific_module{i,3});
	fprintf(fid, '%s\n',Net1_specific_module{i,4});
end
fclose(fid);
filename='Sample2_specific_module.txt';
fid=fopen(filename,'wt');
	fprintf(fid, '%s\t','TF');
	fprintf(fid, '%s\t','TG');
	fprintf(fid, '%s\t','Score');
	fprintf(fid, '%s\n','REs');
for i=1:size(Net2_specific_module,1)
	fprintf(fid, '%s\t',Net2_specific_module{i,1});
	fprintf(fid, '%s\t',Net2_specific_module{i,2});
	fprintf(fid, '%g\t',Net2_specific_module{i,3});
	fprintf(fid, '%s\n',Net2_specific_module{i,4});
end
fclose(fid);
