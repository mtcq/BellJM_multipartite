%Author: Marco Túlio Quintino, https://github.com/mtcq, mtcq.mm@gmail.com
%Requires: nothing
%Last update: 18/Mar/2024

% This script shows that the critical visibility for a pair of Pauli
% measurements to become JM3 is around eta=0.7938;

clear all
X=[0 1;1 0]; Z=[1 0;0 -1]; %Declare Pauli matrices

eta=0.7938;
M(:,:,1)=eta*Z;       
M(:,:,2)=eta*X;
out = IsJM3PairOfMeausrements(M,eta);

eta=0.7937;
M(:,:,1)=eta*Z;       
M(:,:,2)=eta*X;
out = IsJM3PairOfMeausrements(M,eta);

function out = IsJM3PairOfMeausrements(M,eta)

%We start by loading the Sliwa inequalities obtained via faacets.
%The data was obtained from https://github.com/denisrosset/faacets-data/tree/master/solved/L22_22_22.
%The data contenf from faacets is available in the folder L22_22_22 of this repository
%In order to list all inequalities in a systemathic way, I have used the bash script "faacets_converter.sh".
%This script is available in the folder L22_22_22 of this repository

S(:,1) = 	[0, 0, -1, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1,	    0, 0, 0, -1, 0, -1]	;
S(:,2) = 	[0, -1, 0, 1, 0, 1, 0, -1, -1, -1, 0, -1, 0, -1, 1, -1, -1, 0, 0,	    1, 1, -1, -1, 0, 1, 0, -1]	;
S(:,3) = 	[0, 1, 1, 1, -1, 0, -1, 0, -1, 0, -1, -1, -1, 0, 1, 1, -1, -2, 0,	    0, 0, 0, 1, -1, 0, 1, -1]	;
S(:,4) = 	[0, 0, 0, 0, 0, 0, -2, 0, 0, 0, -1, -1, 0, 0, 0, 0, -1, -1, 0, 1,	    -1, 0, 0, 0, 0, 1, -1]	;
S(:,5) = 	[0, 0, 1, 1, -1, 0, 0, -1, -1, -1, 0, -1, 0, -1, 1, -1, -1, 0, 0,	    0, 0, -1, 0, 1, 1, 0, -1]	;
S(:,6) = 	[0, 1, 1, 0, 0, 0, 0, -1, -1, 0, -1, -1, 2, -1, -1, -2, 0, 0, 0, 0,	    0, 0, 1, -1, 0, 1, -1]	;
S(:,7) = 	[0, 0, 0, 0, 1, 1, 0, -1, -1, 0, -2, -2, 2, 0, 0, -2, 0, 0, 0, 0,	    0, 0, 1, -1, 0, 1, -1]	;
S(:,8) = 	[0, 0, 0, 1, 1, 0, -1, -1, 0, 0, -2, 0, 1, 0, -1, -1, 0, -1, 0, 0,	    -2, 0, 1, 1, 0, 1, -1]	;
S(:,9) = 	[0, 0, 0, 0, 1, 1, 0, -1, -1, 0, 0, 0, 0, -2, 0, 0, 0, -2, 0, 0, 0,	    0, 1, -1, 0, 1, -1]	;
S(:,10) = 	[0, 0, 0, 1, 0, 1, -1, 0, -1, 0, 0, 0, 1, -1, -2, -1, -1, 0, 0, 0,	    0, 0, 1, -1, 0, 1, -1]	;
S(:,11) = 	[0, 1, 0, 0, 1, 1, -1, 0, -1, 0, -1, -1, 0, -1, 1, 0, -2, 0, -1, 0,	    -1, 0, 0, 2, -1, 2, -1]	;
S(:,12) = 	[0, 1, 0, 0, 1, 1, -1, 0, -1, 0, -1, -1, 0, 1, -1, 0, 0, -2, -1, 0,	    -1, 0, 0, 2, -1, 2, -1]	;
S(:,13) = 	[0, 1, 0, 0, 1, 1, -1, 0, -1, -1, 0, -1, -1, 0, 1, 0, -2, 0, -2, 1,	    -1, -1, 1, 2, -1, 2, -1]	;
S(:,14) = 	[0, 0, 0, 0, 1, 1, 0, -1, -1, 0, -1, -1, 0, -1, 1, -2, 0, 0, 0, 1,	    1, -2, 0, 0, 0, 1, -1]	;
S(:,15) = 	[0, 2, 0, 0, 0, 0, 0, -2, -2, 0, -1, -1, -1, 2, -1, -1, 1, -2, 0,	    -1, -1, 1, -2, 1, -1, 1, -2]	;
S(:,16) = 	[0, -1, 1, 1, 1, 0, -1, 0, -1, -1, -1, 0, 1, 0, -1, 0, -1, -1, -1,	    0, -1, 0, 1, 1, 1, 1, -2]	;
S(:,17) = 	[0, 0, 0, 2, -1, 1, -2, -1, -1, -2, 0, -2, 1, -2, 1, -1, -2, -1, 0,	    0, -2, -1, 1, 2, 1, 1, -2]	;
S(:,18) = 	[0, 0, 0, 0, 2, 0, 0, 0, -2, 0, 0, 0, 0, -1, -1, 0, -1, -1, 0, 0,	    0, 0, -1, 1, 0, 1, -1]	;
S(:,19) = 	[0, 0, 0, 0, 0, 2, 0, 0, -2, 0, 0, 0, 0, -1, -1, 0, -1, -1, 0, 0,	    0, 0, 1, -1, 0, 1, -1]	;
S(:,20) = 	[0, 0, 0, 0, 1, 1, 0, -1, -1, 0, 0, 0, 0, 0, -2, 0, 0, -2, 0, 0, 0,	    0, -1, 1, 0, 1, -1]	;
S(:,21) = 	[0, 2, 0, 0, 1, 1, 0, -1, -1, 0, -1, -1, 1, -2, 1, 1, -3, 0, 0, -1,	    -1, 1, -1, 2, -1, 2, -1]	;
S(:,22) = 	[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, -3, -1, 0, 0, 0,	    0, 1, -1, 0, 1, -1]	;
S(:,23) = 	[0, 0, 0, 0, 1, 1, 0, -1, -1, 0, 0, 0, -2, 0, 0, -2, 0, 0, 0, 0, 0,	    0, -1, 1, 0, 1, -1]	;
S(:,24) = 	[0, 3, -1, 0, 1, 1, 0, -2, -2, 0, -1, -1, -2, 3, -1, -2, 2, -2, 0,	    -2, -2, 2, -2, 2, -2, 2, -2]	;
S(:,25) = 	[0, -1, -1, 1, 2, 1, -1, -1, -2, 0, 0, 0, -1, -2, -1, -1, -2, -1,	    0, 1, -1, 0, -2, 2, 0, 1, -1]	;
S(:,26) = 	[0, 1, 0, 1, 1, 0, 0, 0, -2, -1, -1, 0, -1, 1, 0, 0, 0, -2, 0, 0,	    -2, 0, 0, 2, -2, 2, 0]	;
S(:,27) = 	[0, 0, 0, -1, 1, 2, -1, -1, -2, 0, -1, -1, 1, -1, -2, -1, -2, -1,	    0, 1, -1, 0, 2, -2, 0, 1, -1]	;
S(:,28) = 	[0, 2, 0, 0, 1, 1, -2, 1, -1, 0, -1, -1, -1, 2, -1, 1, -1, -2, -2,	    1, -1, -1, 1, 2, -1, 2, -1]	;
S(:,29) = 	[0, 0, 0, 0, 1, 1, 0, -1, -1, 0, -1, 1, -1, -1, 0, -1, 0, -1, 0, 1,	    -1, -1, 0, 1, -1, 1, 0]	;
S(:,30) = 	[0, 0, -1, 0, 0, 0, 1, 0, -1, 0, 1, 1, 1, -2, 1, -1, -1, -4, -1, -1,	    0, -1, 2, -1, 0, 1, -3]	;
S(:,31) = 	[0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 1, -1, 0, 0, 0, -2, 1, 1, 0,	    -1, 1, 0, 0, 0, -2]	;
S(:,32) = 	[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0, 0, 0,	    0, 1, 0, 0, -1, 0]	;
S(:,33) = 	[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0,	    0, 0, 1, 0, -1, 0]	;
S(:,34) = 	[0, 0, 0, 1, 0, -1, -1, 0, -1, -1, -1, 0, 0, 2, 0, -1, -1, -2, -1,	    1, 0, 1, -2, -1, 0, 1, -3]	;
S(:,35) = 	[0, 0, 0, 0, 1, -1, -2, -1, -1, -2, -1, -1, -1, 3, 0, -1, -2, -3,	    0, 1, -1, 1, -2, -1, 1, 1, -4]	;
S(:,36) = 	[0, 1, -1, 1, -1, 0, -1, 0, -1, -1, 0, -1, 1, -1, 2, 0, -3, -1, -1,	    1, 0, 0, -2, 0, -1, 1, -2]	;
S(:,37) = 	[0, 0, 0, 0, -1, 1, 0, -1, -1, -1, 0, 1, -1, -1, 2, 0, -3, -1, -1,	    0, -1, 1, -2, 1, 0, 0, -2]	;
S(:,38) = 	[0, 0, 0, 2, -1, -1, 0, -1, -1, 0, -1, -1, -1, 1, 2, -1, -2, -1, 0,	    1, -1, -1, -2, 1, -1, 1, -2]	;
S(:,39) = 	[0, 0, 0, 1, 0, -1, -1, 0, -1, 0, -1, -1, -1, 1, 1, -1, -1, -1, 0,	    1, -1, -1, -1, 1, -1, 1, -1]	;
S(:,40) = 	[0, -1, -1, 1, 1, 0, -1, 0, -1, 0, 0, 0, -2, 0, 0, -2, -2, -2, 0,	    1, -1, -1, -1, 2, -1, 2, -1]	;
S(:,41) = 	[0, 1, -1, -1, 1, 0, -1, 0, -1, -1, 0, -1, -2, 2, -2, -1, -2, -1,	    -1, 1, 0, -1, 1, 2, -2, 0, 0]	;
S(:,42) = 	[0, 0, 0, -1, 0, 1, -1, 0, -1, -1, -1, 0, -1, -1, 1, -1, -1, -1, -1,	    1, 0, -1, 1, 1, -1, 1, -1]	;
S(:,43) = 	[0, -1, -1, 0, 0, 0, 0, -1, 1, 0, -1, 1, -1, -2, -1, -1, -1, -2, 0,	    0, 0, 1, 2, 1, -1, 0, -3]	;
S(:,44) = 	[0, 1, 1, 0, 0, 2, 0, 1, -1, 1, 1, 0, 0, -1, -3, 1, -2, -1, -1, 0,	    1, -2, 3, 1, 1, 1, -4]	;
S(:,45) = 	[0, 0, 0, 1, 2, -1, -3, 2, 1, 3, -1, -2, -1, -4, 1, 2, -3, -3, -1,	    1, -2, 2, 2, 2, -1, 1, -4]	;
S(:,46) = 	[0, 0, 0, 2, -1, 1, -2, -1, 1, 0, -2, 2, -2, 2, -2, -2, -2, 2, 0,	    -2, -2, 0, 3, 1, 0, -1, -3]	;

for i=1:length(S)
    F_Sliwa(:,:,:,i)=faacets2FC_tripartite(S(:,i));
    L_Sliwa(i)=LocalBoundFCtripartite(F_Sliwa(:,:,:,i));
    F_Sliwa_Normalised(:,:,:,i)=F_Sliwa(:,:,:,i);
    F_Sliwa_Normalised(1,1,1,i)=-L_Sliwa(i);
    L_Sliwa_Normalised(i)=LocalBoundFCtripartite(F_Sliwa_Normalised(:,:,:,i));
end

A(:,:,1)=M(:,:,1);    A(:,:,2)=M(:,:,2);
B=A;            C=A;
for i=1:46
    Q(i)=QuantumValueTripartite(F_Sliwa_Normalised(:,:,:,i),A,B,C);
end

[maxQ, ineqMax] = max(Q);
if maxQ>0
    out=1;
    display(['For eta equals to ', num2str(eta)]);
    display(['The these measurements violate the Sliwa inequality number ', num2str(ineqMax)]);
else
    out=0;
    display(['For eta equals to ', num2str(eta)]);
    display(['These measurements do not violate Sliwa inequalities']);
end


end

function ABC=faacets2FC_tripartite(v)

%1 stands for identity
%2 stands for first input
%3 stands for second input

nI=2; %number of inputs
nP=3; %number of parties

ABC=zeros(nI+1,nI+1,nI+1);
for i=1:length(v)
    ABC(i)=v(i);
end
end

function [L_upper L_lower L]=LocalBoundFCtripartite(F)

count=0;
for A1=[-1 1]
    for A2=[-1 1]
        for B1=[-1 1]
            for B2=[-1 1]
                for C1=[-1 1]
                    for C2=[-1 1]
                        count=count+1;
                        A(1,1,1)=A1; A(1,1,2)=A2;
                        B(1,1,1)=B1; B(1,1,2)=B2;
                        C(1,1,1)=C1; C(1,1,2)=C2;
                        L(count)=MakeBellOperatorTripartite(F,A,B,C);
                    end
                end
            end
        end
    end
end
L_upper=max(L);
L_lower=min(L);
end

function BellOp=MakeBellOperatorTripartite(F,A,B,C)

%1 stands for identity
%2 stands for first input
%3 stands for second input

% count the number of inputs
iA=size(A,3);
iB=size(B,3);
iC=size(C,3);
d=size(A,1);

%We now shift the operators to set the very first one as identity
A(:,:,[2 iA+1])=A;
B(:,:,[2 iB+1])=B;
C(:,:,[2 iC+1])=C;

A(:,:,1)=eye(d);
B(:,:,1)=eye(d);
C(:,:,1)=eye(d);
BellOp=0;
for x=1:iA+1
    for y=1:iB+1
        for z=1:iC+1
            BellOp=BellOp+F(x,y,z)*kron(kron(A(:,:,x),B(:,:,y)),C(:,:,z));
        end
    end
end

end

function Q=QuantumValueTripartite(F,A,B,C)
BellOp=MakeBellOperatorTripartite(F,A,B,C);
Q=max(eig(BellOp));
end