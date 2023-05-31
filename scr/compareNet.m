Group1=importdata('Group1');
Group2=importdata('Group2');
if isfile('DesignMat')
     covariates=1;
     Cov_load=importdata('DesignMat');
     [~,sample_idx]=ismember([Group1;Group2],Cov_load.textdata(2:end,1));
     Cov=[ones(length(sample_idx),1) Cov_load.data(sample_idx,:)];
else
     covariates=0;
end
%%%load Exp data and DEG
load('../AllSample/Exp_Opn.mat')
[d,f1]=ismember(Group1,Sample);
[d,f2]=ismember(Group2,Sample);
Exp=Exp(:,[f1;f2]);
Opn=Opn(:,[f1;f2]);
if covariates == 1
    q2_all=dif_test(Exp(:,1:size(Group1,1)),Exp(:,1+size(Group1,1):end),Cov);
else
    q2_all=dif_test(Exp(:,1:size(Group1,1)),Exp(:,1+size(Group1,1):end));
end
q2_all=mafdr(q2_all,'BHFDR',true);
geneName=Symbol(q2_all<0.1);
Exp_all=Exp;
Exp=Exp(q2_all<0.1,:);
q2=q2_all(q2_all<0.1);
[~,filter_idx]=ismember(Symbol,geneName);
%%%%load TRS data
load('../AllSample/TF_binding.mat')
a=sum(TF_binding')';
a1=a(a>0);
TF_binding=TF_binding./(a+1)*median(a1);
load('../AllSample/TFTG_corr.mat')
[d f]=ismember(TFName,Symbol);
TF_max=max(Exp_all(f,:)')';
TF_filter=TF_max>median(TF_max);
R2=R2(TF_filter,q2_all<0.1);
TF_binding=TF_binding(TF_filter,:);
TFName=TFName(TF_filter);
TFExp_median=TFExp_median(TF_filter);
d0=500000;
RE_TG(:,1)=filter_idx(RE_TG(:,1));
RE_TG=RE_TG(RE_TG(:,1)>0,:);
c=double(exp(-1*RE_TG(:,4)/d0).*max(0.2,RE_TG(:,3)));
H1=sparse(RE_TG(:,1),RE_TG(:,2),c,length(geneName),length(Element_name));
for i=1:size(Opn,2)
    BOH(i,:,:)=full(TF_binding.*repmat(Opn(:,i)',size(TF_binding,1),1)*H1');
end
Net1_mean=reshape(mean(BOH(1:size(Group1,1),:,:),1),length(TFName),length(geneName));
Net2_mean=reshape(mean(BOH(1+size(Group1,1):end,:,:),1),length(TFName),length(geneName));
NetPool_reshape=reshape(BOH,size(Group1,1)+size(Group2,1),length(TFName)*length(geneName))';
Net1_mean_fold=Net1_mean./(max(Net1_mean)+0.0001);
Net2_mean_fold=Net2_mean./(max(Net2_mean)+0.0001);
%%%Edge test
if covariates == 1
    [p1,~,NetPool_norm]=dif_test(NetPool_reshape(:,1:size(Group1,1)),NetPool_reshape(:,1+size(Group1,1):end),Cov);
else
    p1=dif_test(NetPool_reshape(:,1:size(Group1,1)),NetPool_reshape(:,1+size(Group1,1):end));
end
p1=mafdr(p1,'BHFDR',true);
p1=max(10^(-16),p1);
p_reshape=-log10(reshape(p1,length(TFName),length(geneName)));
Net1=p_reshape.*(Net1_mean+0.01)./(Net2_mean+0.01).*Net1_mean_fold.*(Net1_mean>Net2_mean);
Net2=p_reshape.*(Net2_mean+0.01)./(Net1_mean+0.01).*Net2_mean_fold.*(Net2_mean>Net1_mean);
%%%%%%%%Figures and Tables
FC=reshape(log2(Net1_mean+0.001./Net2_mean+0.001),length(TFName)*length(geneName),1);
[Net_idx(:,1), Net_idx(:,2)] =ind2sub([length(TFName),length(geneName)],1:length(p1));
diff_idx=(p1<0.1).*(abs(FC)>log2(1.1))==1;
Edge_diff=NetPool_reshape(diff_idx,:);
Edge_name=[TFName(Net_idx(diff_idx,1)) geneName(Net_idx(diff_idx,2))]; 
filename='Diff_TF_TG_pairs.txt';
fid=fopen(filename,'wt');
for i=1:size(Edge_name,1)
	fprintf(fid, '%s\t',Edge_name{i,1});
	fprintf(fid, '%s\n',Edge_name{i,2});
end
fclose(fid);
dlmwrite('Diff_net_raw.txt',[p1(diff_idx) FC(diff_idx) Edge_diff],'\t')
if covariates == 1
    Edge_diff=NetPool_norm(diff_idx,:);
    dlmwrite('Diff_net_normalized.txt',[p1(diff_idx) FC(diff_idx) Edge_diff],'\t')
end
%edgescore=-log10(p1(diff_idx)).*abs(FC(diff_idx));
%[~,f]=sort(edgescore,'descend');
%cg = clustergram(Edge_diff(f(1:min(100,length(f))),:),'Cluster', 'Column','RowPDist', 'correlation', 'Standardize','Row','Colormap',redbluecmap,'ColumnLabels', [Group1;Group2]);
%saveas(gcf,'TF_TG_score_clustering_top100','pdf')
%%%%%%%%%%NodeTest,TG p_value is q2.
[~,f]=ismember(TFName,Symbol);
ExpPool1=Exp_all(f,:);
q1=q2_all(f);
q1(q1<10^(-16))=10^(-16);
q2(q2<10^(-16))=10^(-16);
q1=-log10(q1);
q2=-log10(q2);
%%%Net1 specific 
D1=(Net1>prctile(Net1(Net1>0),90)).*repmat(mean(ExpPool1(:,1:size(Group1,1)),2)>mean(ExpPool1(:,1+size(Group1,1):end),2),1,size(Net1,2));
[a b]=find(D1>0);
c=Net1(D1>0).*q1(a).*q2(b).*abs(R2(D1>0));
Net1_specific=[TFName(a) geneName(b) num2cell(R2(D1>0)) num2cell(q1(a)) num2cell(q2(b)) num2cell(p_reshape(D1>0)) num2cell((Net1_mean(D1>0)+0.01)./(Net2_mean(D1>0)+0.01)) num2cell(Net1_mean_fold(D1>0))];
net_filter=max(q1(a),p_reshape(D1>0))>1;
d=R2(D1>0).*sign(mean(Exp(b,1:size(Group1,1)),2)-mean(Exp(b,1+size(Group1,1):end),2))>0;
Net1_specific=Net1_specific(d.*net_filter>0,:);
[d f]=sort(c(d.*net_filter>0),'descend');
Net1_specific=[Net1_specific(f,:) num2cell(d)];
Net1_specific=Net1_specific(abs(cell2mat(Net1_specific(:,3)))>0.1,:);
w1=double(q1').*(mean(ExpPool1(:,1:size(Group1,1)),2)>mean(ExpPool1(:,1+size(Group1,1):end),2))';
w2=double(q2');
[XNext, YNext,Obj]= FindSubNetwork_bi2(1,0.01,double(Net1),w1,w2);
d1=ismember(Net1_specific(:,1:2),[TFName(XNext>prctile(XNext,90));geneName(YNext>prctile(YNext,90))]);
Net1_specific_module=Net1_specific(d1(:,1).*d1(:,2)>0,:);
%%Net2 specific 
D1=(Net2>prctile(Net2(Net2>0),90)).*repmat(mean(ExpPool1(:,1:size(Group1,1)),2)<mean(ExpPool1(:,1+size(Group1,1):end),2),1,size(Net1,2));
[a b]=find(D1>0);
c=Net2(D1>0).*q1(a).*q2(b).*abs(R2(D1>0));
Net2_specific=[TFName(a) geneName(b) num2cell(R2(D1>0)) num2cell(q1(a)) num2cell(q2(b)) num2cell(p_reshape(D1>0)) num2cell((Net2_mean(D1>0)+1)./(Net1_mean(D1>0)+1)) num2cell(Net2_mean_fold(D1>0))];
net_filter=max(q1(a),p_reshape(D1>0))>1;
d=R2(D1>0).*sign(mean(Exp(b,1+size(Group1,1):end),2)-mean(Exp(b,1:size(Group1,1)),2))>0;
Net2_specific=Net2_specific(d.*net_filter>0,:);
[d f]=sort(c(d.*net_filter>0),'descend');
Net2_specific=[Net2_specific(f,:) num2cell(d)];
Net2_specific=Net2_specific(abs(cell2mat(Net2_specific(:,3)))>0.1,:);
w1=double(q1').*(mean(ExpPool1(:,1:size(Group1,1)),2)'<mean(ExpPool1(:,1+size(Group1,1):end),2)');
w2=double(q2');
[XNext, YNext,Obj]= FindSubNetwork_bi2(1,0.01,double(Net2),w1,w2);
d1=ismember(Net2_specific(:,1:2),[TFName(XNext>prctile(XNext,90));geneName(YNext>prctile(YNext,90))]);
Net2_specific_module=Net2_specific(d1(:,1).*d1(:,2)>0,:);
%%%%%%%%%%%%%%%%common
Net3=sqrt(Net1_mean/max(max(Net1_mean)).*Net2_mean/max(max(Net2_mean)));
D1=Net3>max(prctile(Net3(Net3>0),95),0.01);
[a b]=find(D1>0);
c=Net1_mean(D1>0)/max(max(Net1_mean)).*Net2_mean(D1>0)/max(max(Net2_mean));
Net_common=[TFName(a) geneName(b) num2cell(R2(D1>0)) num2cell(Net1_mean(D1>0)/max(max(Net1_mean))) num2cell(Net2_mean(D1>0)/max(max(Net2_mean)))];
[d f]=sort(c,'descend');
Net_common=[Net_common(f,:) num2cell(d)];
Net_common=Net_common(abs(cell2mat(Net_common(:,3)))>0.1,:);
[XNext, YNext,Obj]= FindSubNetwork_bi1(2,double(Net3));
d1=ismember(Net_common(:,1:2),[TFName(XNext>prctile(XNext,95));geneName(YNext>prctile(YNext,95))]);
Net_common_module=Net_common(d1(:,1).*d1(:,2)>0,:);
%%%%%%%%%%%%%%%%write
Node_label=1*(mean(Exp(:,1:size(Group1,1)),2)>mean(Exp(:,1+size(Group1,1):end),2));
Node_label(Node_label==0)=2;
filename='Node_label.txt';
fid=fopen(filename,'wt');
for i=1:size(geneName,1)
	fprintf(fid, '%s\t',geneName{i,1});
	fprintf(fid, '%g\n',Node_label(i,1));
end
fclose(fid);
filename='Group1_specific_network.txt';
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
filename='Group2_specific_network.txt';
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
filename='Group1_specific_module.txt';
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
filename='Group2_specific_module.txt';
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
filename='Group1_Group2_common_network.txt';
fid=fopen(filename,'wt');
	fprintf(fid, '%s\t','TF');
	fprintf(fid, '%s\t','TG');
	fprintf(fid, '%s\t','corr');
	fprintf(fid, '%s\t','Activity_Group1');
	fprintf(fid, '%s\t','Activity_Group2');
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
filename='Group1_Group2_common_module.txt';
fid=fopen(filename,'wt');
	fprintf(fid, '%s\t','TF');
	fprintf(fid, '%s\t','TG');
	fprintf(fid, '%s\t','corr');
	fprintf(fid, '%s\t','Activity_Group1');
	fprintf(fid, '%s\t','Activity_Group2');
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

filename='Group1_specific_TFnetwork.txt';
Net1_TFspecific=Net1_specific(ismember(Net1_specific(:,2),TFName),:);
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
filename='Group2_specific_TFnetwork.txt';
Net2_TFspecific=Net2_specific(ismember(Net2_specific(:,2),TFName),:);
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
