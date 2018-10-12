function [XNext, YNext,Obj]= FindSubNetwork_bi2(beta,Lambda,SPMatrix,F1,F2)
    [n1,n2]=size(SPMatrix);
    IterMax=1000;
    Epsilon=0.0001;
    Threshold=1e-8;
    p=50000;
    Iter=1;        %The number of the iterative
    %X0=ones(n,1)*(1/n);
    for i=1:5
    X0=sparse(rand(1,n1));
    Y0=sparse(rand(1,n2));
    alpha= (X0*SPMatrix*Y0'+Lambda*(F1*X0'+F2*Y0'))/beta;       % Langrange mulitiplier
    XPre= ( X0.*(Y0*SPMatrix'+Lambda*F1)/(beta*alpha) ).^(1/beta);
    YPre= ( Y0.*(XPre*SPMatrix+Lambda*F2)/(beta*alpha) ).^(1/beta);
    Error=max(norm(XPre-X0,'fro')+norm(YPre-Y0,'fro'),10*Epsilon);
    YNext=YPre;
    while (Iter < IterMax ) && (Error > Epsilon)  
        alpha= (XPre*SPMatrix*YPre'+Lambda*(F1*XPre'+F2*YPre'))/beta;
        XNext= ( XPre.*(YPre*SPMatrix'+Lambda*F1)/(beta*alpha) ).^(1/beta);
        YNext= ( YPre.*(XPre*SPMatrix+Lambda*F2)/(beta*alpha) ).^(1/beta);
        Error=norm(XNext-XPre,'fro')+norm(YNext-YPre,'fro');
        Iter=Iter+1;
        XPre=XNext;
        YPre=YNext;
        Obj=trace(XPre*SPMatrix*YPre');
    end 
    ObjValue(i,1)=Obj;
    XNextPool{1,i}=XNext;
    YNextPool{1,i}=YNext;
    end
    [d f]=max(ObjValue);
    XNext=XNextPool{1,f};
    YNext=YNextPool{1,f};