function [p1,p2,Res]=dif_test(A,B,Cov)
if nargin > 2
AB=[A,B];
Res= AB-AB*Cov*pinv(Cov'*Cov)*Cov';
A=Res(:,1:size(A,2));
B=Res(:,1+size(A,2):end);
end
%p2 ttest2 pval
ds1=size(A,2)-1;
ds2=size(B,2)-1;
ds=size(A,2)+size(B,2)-1;
n=size(A,1);
s1=std(A')';
s2=std(B')';
s1(s1<0.001)=0.001;
s2(s2<0.001)=0.001;
t=((mean(A')'-mean(B')')./sqrt((s1.^2)/(ds1+1)+(s2.^2)/(ds2+1)));
p2=2*(1-tcdf(abs(t),ds));
z=log([s1;s2].^2);
e=z-psi(ds/2)+log(ds/2);
inv_d0=double(mean((e-mean(e)).^2*n/(n-1)-psi(1,ds/2)));
d0 = 2*fsolve(@(x)psi(1,x)-inv_d0,0.1,statset('display','off'));
ss=exp(mean(e)+psi(d0/2)-log(d0/2));
ss1=sqrt((d0*ss+ds1*s1.^2)/(ds1+d0));
ss2=sqrt((d0*ss+ds2*s2.^2)/(ds2+d0));
t1=((mean(A')'-mean(B')')./sqrt((ss1.^2)/(ds1+1)+(ss2.^2)/(ds2+1)));
p1=2*(1-tcdf(abs(t1),ds));
