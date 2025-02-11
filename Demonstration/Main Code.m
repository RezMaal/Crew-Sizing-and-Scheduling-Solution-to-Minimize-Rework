%% Crew Sizing and Scheduling
%% Initializations
D=readmatrix("Distance.txt");
Arr=readmatrix("Arrangement.txt");
Space_BIM=readmatrix('space_BIM.txt');
Err_BIM=readmatrix('err_BIM.txt');
Score_BIM=zeros(length(Err_BIM),1);
az=unique(Space_BIM(:,1));
for i=1:length(Err_BIM)
    x = Space_BIM(Space_BIM(:,1)==az(i),2); % X-coordinates
    y = Space_BIM(Space_BIM(:,1)==az(i),3); % Y-coordinates
    Score_BIM(i,1)=Err_BIM(Err_BIM(:,1)==az(i),2)^2/polyarea(x,y)^1.5*3;
end
x_min = min(Score_BIM);
x_max = max(Score_BIM);
a=8; b=25;
% Avoid division by zero if all values in x are the same
if x_max == x_min
    Score_BIM_C = a * ones(size(Score_BIM)); % If all values are the same, set all to 'a'
else
    Score_BIM_C = a + ((Score_BIM - x_min) / (x_max - x_min)) * (b - a);
end
Score_BIM_C=ceil(Score_BIM_C);

%% Input Parameters
k=10; n=length(az);
M_col=15*ones(1,15);
[~,indices] = ismember(Arr,az-8);
alpha_col=Score_BIM_C(indices);
M_c=1.5*max(alpha_col)*ones(k,1);
C_c=120; C_o=100; C_t=70;
Q=D-1;

%% Crew Adjustment No Space No Takt
result=Only_Crew(k,n,C_c/max(M_col),C_o,alpha_col,M_col,M_c);
X_C=reshape(result.x(1:n*k),[k,n]);
% estimating time and cost
ss=sum(X_C,2);
[~,a]=sort(ss,'descend');
xx_C=X_C(a,:);
t_C=sum(ss>0);
Co_C=sum(X_C,"All")*C_c/max(M_col)+C_o*t_C;

%% Crew Adjustment with Space No Takt
result=Only_Crew_Spatial(k,n,C_c/max(M_col),C_o,alpha_col,M_col,M_c,Q);
X_Cs=reshape(result.x(1:n*k),[k,n]);
% estimating time and cost
ss=sum(X_Cs,2);
[~,a]=sort(ss,'descend');
xx_Cs=X_Cs(a,:);
t_Cs=sum(ss>0);
Co_Cs=sum(X_Cs,"All")*C_c/max(M_col)+C_o*t_Cs;

%% Crew Adjustment No Space with Takt
result=Takt_Crew(k,n,C_c/max(M_col),C_o,C_t,alpha_col,M_col,M_c);
X_Ct=reshape(result.x(1:n*k),[k,n]);
% estimating time and cost
ss=sum(X_Ct,2);
[~,a]=sort(ss,'descend');
xx_Ct=X_Ct(a,:);
t_Ct=sum(ss>0);
Co_Ct=sum(X_Ct,"All")*C_c/max(M_col)+C_o*t_Ct;

%% Crew Adjustment with Space with Takt
result=Takt_Crew_Spatial(k,n,C_c/max(M_col),C_o,C_t,alpha_col,M_col,M_c,Q);
X_Cts=reshape(result.x(1:n*k),[k,n]);
% estimating time and cost
ss=sum(X_Cts,2);
[~,a]=sort(ss,'descend');
xx_Cts=X_Cts(a,:);
t_Cts=sum(ss>0);
Co_Cts=sum(X_Cts,"All")*C_c/max(M_col)+C_o*t_Cts;

MX=max(max([xx_C;xx_Ct;xx_Cs;xx_Cts])); mX=min(min([xx_C;xx_Ct;xx_Cs;xx_Cts])); 
