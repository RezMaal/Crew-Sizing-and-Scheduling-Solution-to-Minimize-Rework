function [Attribute,Bounds,Area,Mean]=readSTLPolyGraph(STL)
    TR = stlread(STL);
    N = faceNormal(TR);
    P=TR.Points;
    Con=TR.ConnectivityList;
    Cn=size(Con,1);
    for i=1:Cn
        p=P(Con(i,:),:);
        mp=mean(p);
        N(i,4)=-dot(mp,N(i,1:3));
    end
    Dn=round(N(:,1:3)*N(:,1:3)',10);   
    Adj=zeros(Cn,Cn);
    for i=1:Cn
        aa=ismember(Con(i+1:end,:),Con(i,:));
        aa=sum(aa,2);
        Adj(i,i+find(aa==2))=1;
    end
    Adj=Adj+Adj';
    Adj(abs(Dn)~=1)=0;
    G=graph(Adj);
    bins = conncomp(G);
    mx=max(bins);
    Attribute=zeros(mx,4);
    Bounds=cell(mx,1);
    Mean=zeros(mx,3);
    Area=zeros(mx,1);
    warning('off','all')
    for i=1:mx
        ii=find(bins==i);
        Rt=Rodrigues(N(ii(1),1:3),[0,0,1]);
        In=unique(Con(ii,:));
        P_Rot=P(In,:)*Rt;
        Index=boundary(P_Rot(:,1:2),0);
        Pri=P_Rot(Index,1:2);
        pgon=polyshape(Pri);
        Pr=pgon.Vertices;
        inx=knnsearch(Pri,Pr);
        PB=P(In(Index(inx)),:);
        Bounds{i,1}=PB;
        [~,d]=max(abs(N(ii(1),1:3)));
        if sign(N(ii(1),d))<0
            V=-N(ii(1),:);
        else
            V=N(ii(1),:);
        end
        Attribute(i,:)=V;
        Mean(i,:)=mean(PB);
        Area(i,1)=area(pgon);
    end
end
% Rodrigues's formula for Vector Rotation: This function provides the
% rotation matrix to rotate vector N to V. N and V are 1x3 matrices
function Rot=Rodrigues(N,V)
    teta=acos(dot(N,V)/norm(N)/norm(V));
    k=cross(N,V);
    if norm(k)==0
        Rot=eye(length(N),length(N));
    else
        k=k/norm(k);
        K=[0,-k(3),k(2);k(3),0,-k(1);-k(2),k(1),0];
        R=eye(3)+sin(teta)*K+(1-cos(teta))*K^2;
        Rot=R';
    end
end