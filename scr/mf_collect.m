filename='knownResults.txt';
fileID = fopen(filename);
C = textscan(fileID,'%s %s %f %f %f %f %f %s %f %f %s','HeaderLines', 1);
fclose(fileID);
Score=[-log10(C{1,3}) (C{1,7}+0.1)./(C{1,10}+0.1)];
Score(Score(:,1)>100,1)=100;
Score(:,3)=sqrt(Score(:,1).*Score(:,2));
[d f]=sort(Score(:,3),'descend');
Score=Score(f,:);
Name=C{1,1}(f);
filename='knownResults_rank.txt';
fid=fopen(filename,'wt');
	fprintf(fid, '%s\t','Motif');
	fprintf(fid, '%s\t','-log10(p)');
	fprintf(fid, '%s\t','FoldChange');
	fprintf(fid, '%s\n','Score');
for i=1:size(Name,1)
	fprintf(fid, '%s\t',Name{i,1});
	fprintf(fid, '%g\t',Score(i,1));
	fprintf(fid, '%g\t',Score(i,2));
	fprintf(fid, '%g\n',Score(i,3));
 end
 fclose(fid);
%%%%%%%%%%%%%%%
load('../../../Data/MotifMatch_species_rmdup.mat')
[d f]=ismember(Match2(:,1),Name);
Score1=Score(f(d==1),:);
Match2=Match2(d==1,:);
[TF ia ic]=unique(Match2(:,2));
TFScore=accumarray(ic,Score1(:,3),[],@max);
[d f]=sort(TFScore,'descend');
Results=[TF(f) num2cell(TFScore(f))];
filename='knownResults_TFrank.txt';
fid=fopen(filename,'wt');
for i=1:size(Results,1)
	fprintf(fid, '%s\t',Results{i,1});
	fprintf(fid, '%g\n',Results{i,2});
 end
 fclose(fid);