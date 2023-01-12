fileID = fopen('openness2.bed');
C = textscan(fileID,'%s %f32 %f32');
fclose(fileID);
Element_name=C{1,1};
Opn=double(C{1,2});
Opn_median=double(C{1,3});
load('../../Data/MotifMatch_mouse_rmdup.mat')
load('TFTG_corr.mat')
TF_binding=mfbs(TFName,Element_name,motifName,motifWeight,Match2);
fileID = fopen('../../Input/toreplace.txt');
C = textscan(fileID,'%s %f32');
fclose(fileID);
Symbol=C{1,1};
G=C{1,2};
%%%%CR binding
load('../../Data/CRInfo_mouse.mat')
eita0=-30.4395;
eita1=0.8759;
[d f]=ismember(C_TFName,TFName);
C_TFName=C_TFName(d==1,:);
TFS=TFS(d==1,:);
CR_TF=CR_TF(:,d==1);
TFB=full(TF_binding(f(d==1),:));
C_TFExp=zeros(length(C_TFName),1);
[d f]=ismember(C_TFName,Symbol);
C_TFExp(d==1,:)=log2(1+G(f(d==1),:));
[d1 f1]=ismember(C_TFName,TFName);
TFBO=(repmat(C_TFExp,1,length(Opn)).*repmat(C_TFExp./TFS,1,length(Opn)).*TF_binding(d1==1,:).*repmat(Opn',length(C_TFName),1)).^0.25;
CRB=eita0+eita1*CR_TF*TFBO;
CRB_P=1-exp(CRB')./(1+exp(CRB'));
filename='CR_binding_pval.txt';
 fid=fopen(filename,'wt');
for i=1:length(CRName)-1
         fprintf(fid, '%s\t',CRName{i,1});
end
fprintf(fid, '%s\n',CRName{i+1,1});
 fclose(fid);
dlmwrite('CR_binding_pval.txt',CRB_P,'-append','delimiter','\t')
%%%%%%%%%%%%
alhfa=0.5;
Opn_median=log2(1+Opn_median);
Opn1=log2(1+Opn);
Opn=double(Opn1.*(Opn1./(Opn_median+0.5)));
geneName=intersect(List,Symbol);
[d f]=ismember(geneName,List);
R2=R2(:,f);
Exp_median=Exp_median(f);
[d f]=ismember(geneName,Symbol);
G=G(f);
[d1 f1]=sort(G);
[d2 f2]=sort(Exp_median);
G1(f1,1)=d2;
G=(G1.^(alhfa)).*(G1./(Exp_median+0.5));
[d f]=ismember(TFName,geneName);
TFName=TFName(d==1,:);
TF_binding=TF_binding(d==1,:);
TFExp=G(f(d==1));
R2=R2(d==1,:);
fileID = fopen('./Enrichment/knownResults_TFrank.txt');
C = textscan(fileID,'%s %f32');
fclose(fileID);
[d f]=ismember(TFName,C{1,1});
TF_motif=zeros(length(TFName),1);
TF_motif(d==1)=C{1,2}(f(d==1));
TFExp=TFExp.*(TF_motif);

fileID = fopen('peak_gene_100k_corr.bed');
C = textscan(fileID,'%s %s %f32 %f32');
fclose(fileID);

[d f]=ismember(C{1,1},Element_name);
[d1 f1]=ismember(C{1,2},geneName);
[f2 ia ic]=unique([f(d.*d1==1) f1(d.*d1==1)],'rows');
c3=accumarray(ic,C{1,3}(d.*d1==1),[],@min);
c4=accumarray(ic,C{1,4}(d.*d1==1),[],@min);
RETG_w=c4;
RETG_name=[Element_name(f2(:,1))  geneName(f2(:,2))];
c4(c4<0.2)=0;
d0=500000;
c=double(exp(-1*c3/d0).*c4);
H1=sparse(f2(:,2),f2(:,1),c,length(geneName),length(Element_name));
TFO=TF_binding.*repmat(Opn',size(TF_binding,1),1);
BOH=TFO*H1';
Score=(TFExp*G').*(2.^abs(R2)).*full(BOH);
Score(isnan(Score))=0;
dlmwrite('TFTG_regulationScore.txt',Score,'\t')
Score_norm=(2.^abs(R2)).*full(BOH);
Score_norm(isnan(Score_norm))=0;
dlmwrite('TFTG_regulationScore_norm.txt',Score_norm,'\t')
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
load('../../Data/TFTG_mouse_nagetriveControl.mat')
[d f]=ismember(Back_net(:,1),TFName);
[d1 f1]=ismember(Back_net(:,2),geneName);
f2=[f(d.*d1==1) f1(d.*d1==1)];
Back_score=Score((f2(:,2)-1)*size(Score,1)+f2(:,1));
Cut=prctile(Back_score,99);
[a b]=find((Score>Cut)==1);
c=find((Score>Cut)==1);
c1=full(Score(c));
Net=[TFName(a) geneName(b)];

d=ismember(RETG_name(:,2),unique(Net));
RETG_name=RETG_name(d,:);
RETG_w=RETG_w(d,:);
filename='toreplace_RE_TG.txt';
fid=fopen(filename,'wt');
	fprintf(fid, '%s\t','RE');
	fprintf(fid, '%s\t','TG');
	fprintf(fid, '%s\n','Weight');
for i=1:size(RETG_w,1)
	fprintf(fid, '%s\t',RETG_name{i,1});
	fprintf(fid, '%s\t',RETG_name{i,2});
	fprintf(fid, '%g\n',RETG_w(i,1));
end
fclose(fid);

[a1,a2]=sort(H1','descend');
a1=a1(1:10,:);
a2=a2(1:10,:);
TFTG_RE=arrayfun(@(i)  strjoin(Element_name(a2((TFO(a(i),a2(:,b(i)))>0)'.*(a1(:,b(i))>0)==1,b(i)))',';'),[1:length(a)]','UniformOutput',false);
[d f]=sort(c1,'descend');
Net=[Net(f,:) num2cell(d) TFTG_RE(f)];
filename='toreplace_network.txt';
fid=fopen(filename,'wt');
	fprintf(fid, '%s\t','TF');
	fprintf(fid, '%s\t','TG');
	fprintf(fid, '%s\t','Score');
	fprintf(fid, '%s\t','FDR');
	fprintf(fid, '%s\n','REs');
for i=1:size(Net,1)
	fprintf(fid, '%s\t',Net{i,1});
	fprintf(fid, '%s\t',Net{i,2});
	fprintf(fid, '%g\t',Net{i,3});
	fprintf(fid, '%g\t',(sum(Back_score>d(i))+1)/length(Back_score));
	fprintf(fid, '%s\n',Net{i,4});
end
fclose(fid);
