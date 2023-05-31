%%%%%%%load data
Sample=importdata("SampleNameFile")
for i=1:size(Sample,1)
    fileID = fopen(['../../Input/',Sample{i,1},'.txt']);
    C = textscan(fileID,'%s %f32');
    fclose(fileID);
    Exp(:,i)=double(C{1,2});
    Symbol=C{1,1};
    fileID = fopen(['openness_',Sample{i,1},'.bed']);
    C = textscan(fileID,'%s %f32');
    fclose(fileID);
    Opn(:,i)=double(C{1,2});
    PeakName=C{1,1};
end
Exp=log2(1+quantilenorm(Exp));
Opn=log2(1+quantilenorm(Opn));
%%%%%%%%%%%RE_TG
fileID = fopen('peak_gene.txt');
C = textscan(fileID,'%s %s');
fclose(fileID);
[d1 f1]=ismember(C{1,2},Symbol);
A1=zeros(length(d1),size(Exp,2));
A1(d1==1,:)=Exp(f1(d1==1),:);
[d2 f2]=ismember(C{1,1},PeakName);
Opn1=Opn(f2,:);
Z1=zscore(A1');
Z2=zscore(Opn1');
r=sum(Z1.*Z2)/(size(Z1,1)-1);
r=r';
dlmwrite('corr',r,'\t')
dis=dlmread('dis');
RE_TG=[f1(d1) f2(d1) r(d1) dis(d1)];
Element_name=PeakName;
save('Exp_Opn.mat','Exp','Opn','Symbol','Element_name','Sample','RE_TG')
%%%%TF-TG
TFName=importdata('../../Data/TFName_speciesFull.txt');
[d,f]=ismember(TFName,Symbol);
TFExp=Exp(f(d==1),:);
TFName=TFName(d==1);
R2=corr(TFExp',Exp');
List=Symbol;
Exp_median=median(Exp')';
TFExp_median=median(TFExp')';
save('TFTG_corr.mat','TFName','List','R2','Exp_median','TFExp_median')
%%%%%%motif_activity
load('../../Data/MotifMatch_speciesFull_rmdup.mat')
[~,Motif_binding]=mfbs(TFName,PeakName,motifName,motifWeight,Match2);
Opn_median=median(Opn,2);
Motif_activity=Motif_binding*(Opn./(Opn_median+0.5));
[Motif_TF,p]=corr(Motif_activity',TFExp');
p(isnan(Motif_TF))=1;
Motif_TF(isnan(Motif_TF))=0;
[~,f1]=ismember(Match2(:,1),motifName);
[~,f2]=ismember(Match2(:,2),TFName);
d=f1.*f2>0;
Motif_TF_map=sparse(f1(d),f2(d),ones(sum(d),1),length(motifName),length(TFName));
map_id=[f1(d) f2(d)];
map_id(:,3)=Motif_TF((map_id(:,2)-1)*length(motifName)+map_id(:,1));
[d f]=sort(map_id(:,3),'descend');
Map_score=[motifName(map_id(f,1)) TFName(map_id(f,2)) num2cell(map_id(f,3))];
filename='TF_motif_corr.txt';
fid=fopen(filename,'wt');
for i=1:size(Map_score,1)
fprintf(fid, '%s\t',Map_score{i,1});
fprintf(fid, '%s\t',Map_score{i,2});
fprintf(fid, '%g\n',Map_score{i,3});
end
fclose(fid);
Map_matrix=(Motif_TF>0).*(p<0.1).*(Motif_TF_map>0);
Map_matrix_update=Map_matrix;
for i=1:size(Map_matrix,2)
  if(sum(Map_matrix(:,i))==0)
    a=find(Motif_TF_map(:,i)>0);
    [~,idx]=max(Motif_TF(a,i));
    Map_matrix_update(a(idx),i)=1;
  end
end
TF_binding = sparse(Map_matrix_update'*Motif_binding);
save('TF_binding.mat','TFName','TF_binding')
