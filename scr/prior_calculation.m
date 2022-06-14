%%%%%%%load data
Sample=importdata(SampleNameFile)
for i=1:size(Sample,1)
    Exp(:,i)=dlmread(['../../Input/',Sample{i,1},'.txt'],0,1);
    Opn(:,i)=dlmread(['openness_',Sample{i,1},'.bed'],0,1);
end
Exp=log2(1+quantilenorm(Exp));
Opn=log2(1+quantilenorm(Opn));
fileID = fopen(['../../Input/',Sample{1,1},'.txt']);
C = textscan(fileID,'%s %*[^\n]');
fclose(fileID);
Symbol=C{1,1};
fileID = fopen(['openness_',Sample{1,1},'.bed']);
C = textscan(fileID,'%s %*[^\n]');
fclose(fileID);
PeakName=C{1,1};
%%%%%%%%%%%RE_TG
 fileID = fopen('peak_gene.txt');
C = textscan(fileID,'%s %s');
fclose(fileID);
[d f]=ismember(C{1,2},Symbol);
A1=zeros(length(d),size(A,2));
A1(d==1,:)=A(f(d==1),:);
[d f]=ismember(C{1,1},PeakName);
Opn1=Opn(f,:);
Z1=zscore(A1');
Z2=zscore(Opn1');
r=sum(Z1.*Z2)/(size(Z1,1)-1);
r=r';
dlmwrite('corr',r,'\t')
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