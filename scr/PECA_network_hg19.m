fileID = fopen('region.txt');
C = textscan(fileID,'%s %s %s  %s');
fclose(fileID);
Element_name=C{1,4};
load('../../Data/MotifMatch_human_rmdup.mat')
load('../../Data/TFTG_corr_human.mat')
TF_binding=mfbs(TFName,Element_name,motifName,motifWeight,Match2);
Opn=dlmread('openness.bed','\t',0,1);
fileID = fopen('../../Input/toreplace.txt');
C = textscan(fileID,'%s %f32');
fclose(fileID);
Symbol=C{1,1};
G=C{1,2};
geneName=intersect(List,Symbol);
[d f]=ismember(geneName,List);
R2=R2(:,f);
[d f]=ismember(geneName,Symbol);
G=G(f);
[d f]=ismember(TFName,geneName);
TFName=TFName(d==1,:);
TF_binding=TF_binding(d==1,:);
TFExp=G(f(d==1));

fileID = fopen('peak_gene_100k_corr.bed');
C = textscan(fileID,'%s %s %f32 %f32');
fclose(fileID);

[d f]=ismember(C{1,1},Element_name);
[d1 f1]=ismember(C{1,2},geneName);
[f2 ia ic]=unique([f(d.*d1==1) f1(d.*d1==1)],'rows');
c3=accumarray(ic,C{1,3}(d.*d1==1),[],@min);
c4=accumarray(ic,C{1,4}(d.*d1==1),[],@min);
c4(c4<0)=0;
d0=100000;
c=double(exp(-1*c3/d0).*c4);
H1=sparse(f2(:,2),f2(:,1),c,length(geneName),length(Element_name));
TFO=TF_binding.*Opn';
BOH=TFO*H1';
Score=sqrt(TFExp*G').*(2.^abs(R2)).*full(BOH);
dlmwrite('TFTG_regulationScore.txt',Score,'\t')
filename='TFName.txt';
 fid=fopen(filename,'wt');
for i=1:size(TFName,1)
 fprintf(fid, '%s\n',TFName{i,1});
 end
 fclose(fid);

filename='TGName.txt';
 fid=fopen(filename,'wt');
for i=1:size(geneName)
 fprintf(fid, '%s\n',geneName{i,1});
 end
 fclose(fid);

%%%%%%%%%
load('../../Data/TFTG_human_nagetriveControl.mat')
[d f]=ismember(Back_net(:,1),TFName);
[d1 f1]=ismember(Back_net(:,2),geneName);
f2=[f(d.*d1==1) f1(d.*d1==1)];
Back_score=Score((f2(:,2)-1)*size(Score,1)+f2(:,1));
Cut=prctile(Back_score,99.9);
[a b]=find(Score>Cut);
c=find(Score>Cut);
c1=full(Score(c));
Net=[TFName(a) List(b)];
TFTG_RE=arrayfun(@(i)  strjoin(Element_name(find((TFO(a(i),:)>0).*(H1(b(i),:)>0)==1))',';'),[1:length(a)]','UniformOutput',false);
[d f]=sort(c1,'descend');
Net=[Net(f,:) num2cell(d) TFTG_RE(f)];
filename='toreplace_network.txt';
fid=fopen(filename,'wt');
	fprintf(fid, '%s\t','TF');
	fprintf(fid, '%s\t','TG');
	fprintf(fid, '%s\t','Score');
	fprintf(fid, '%s\n','REs');
for i=1:size(Net,1)
	fprintf(fid, '%s\t',Net{i,1});
	fprintf(fid, '%s\t',Net{i,2});
	fprintf(fid, '%g\t',Net{i,3});
	fprintf(fid, '%s\n',Net{i,4});
end
fclose(fid);
