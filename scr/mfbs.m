function TF_binding=mfbs(TFName,Element_name,motifName,motifWeight,Match2)
    fileID = fopen('MotifTarget.txt');
    C = textscan(fileID,'%s %s %f32');
    fclose(fileID);
    f3=C{1,3};
    [d1 f1]=ismember(C{1,1},Element_name);
    [d2 f2]=ismember(C{1,2},motifName);
    t1=setdiff([1:length(motifName)],unique(f2));
    f2=[f2((d1.*d2)==1);t1'];
    f1=[f1((d1.*d2)==1);ones(length(t1),1)];
    f3=[f3((d1.*d2)==1);zeros(length(t1),1)];
    t1=setdiff([1:length(Element_name)],unique(f1));
    f1=[f1;t1'];
    f2=[f2;ones(length(t1),1)];
    f3=[f3;zeros(length(t1),1)];
    Motif_binding=sparse(double(f2),double(f1),double(f3),length(motifName),size(Element_name,1));
    Motif_binding=diag(1./(motifWeight+0.1))*Motif_binding;
    Motif_binding=log(1+Motif_binding);

TF_binding=zeros(length(TFName),length(Element_name));
[d f1]=ismember(Match2(:,1),motifName);
[d f2]=ismember(Match2(:,2),TFName);
for i=1:length(TFName)
	a=find(f2==i);
	if length(a)>1
	TF_binding(i,:)=max(Motif_binding(f1(a),:));
	else if length(a)==1
	TF_binding(i,:)=Motif_binding(f1(a),:);
	else
	TF_binding(i,:)=zeros(1,length(Element_name));
	end
	end
end
TF_binding=sparse(TF_binding);