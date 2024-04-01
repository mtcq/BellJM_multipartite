%Author: Marco Túlio Quintino, https://github.com/mtcq, mtcq.mm@gmail.com
%Requires: nothing
%Last update: 01/Apr/2024

% This script shows that the critical visibility when two parties performs the measurements T presented in:
% E. Bene and T. Vértesi, Measurement incompatibility does not give rise to Bell violation in general, New Journal of Physics 20, 013021 (2018), arXiv:1705.10069

%We find that JM_2(T) is precisely greater than eta = sqrt(2/sqrt(7)) \approx 0.8694
% For that, we evaluate the CHSH value analytically, to find the critical  visibility if  sqrt(2/sqrt(7)).
% Then we show that, with this visibility, one cannot violate the I3322 inequality using the T measurements

clear all
Id=eye(2); % X=[0 1;1 0]; Y=[0 -sqrt(-1);sqrt(-1) 0]; Z=[1 0;0 -1]; H=[1 1; 1 -1]/sqrt(2);
% eta=.8944;
eta=.87;

a0=[cos(2*0*pi/3) 0 sin(2*0*pi/3)];
a1=[cos(2*1*pi/3) 0 sin(2*1*pi/3)];
a2=[cos(2*2*pi/3) 0 sin(2*2*pi/3)];
T(:,:,1,1)=Bloch2Rho(a0);
T(:,:,1,2)=Bloch2Rho(a1);
T(:,:,1,3)=Bloch2Rho(a2);
for x=1:3
    T(:,:,2,x)=Id-T(:,:,1,x);
    Teta(:,:,1,x)=eta*T(:,:,1,x)+(1-eta)*Id/2;
    Teta(:,:,2,x)=Id-Teta(:,:,1,x);
end
for x=1:3
    T_ob(:,:,x)=T(:,:,1,x)-T(:,:,2,x);
    T_ob_eta(:,:,x)=eta*T_ob(:,:,x);
end

I3322_FC=[0 1 1 0;-1 1 1 1;-1 1 1 -1;0 1 -1 0];
% I3322_FC(1,1)=-4;


%+++  1
I3322Op=MakeBellOperatorBipartite(I3322_FC,T_ob_eta,T_ob_eta);
MasValueI3322(1)=max(eig(I3322Op));
%++-  2
T_ob_eta(:,:,3)=-T_ob_eta(:,:,3);
I3322Op=MakeBellOperatorBipartite(I3322_FC,T_ob_eta,T_ob_eta);
MasValueI3322(2)=max(eig(I3322Op));
%+-+  3
T_ob_eta(:,:,2)=-T_ob_eta(:,:,2);
T_ob_eta(:,:,3)=-T_ob_eta(:,:,3);
I3322Op=MakeBellOperatorBipartite(I3322_FC,T_ob_eta,T_ob_eta);
MasValueI3322(3)=max(eig(I3322Op));
%+-- 4
T_ob_eta(:,:,2)=-T_ob_eta(:,:,2);
T_ob_eta(:,:,3)=-T_ob_eta(:,:,3);
I3322Op=MakeBellOperatorBipartite(I3322_FC,T_ob_eta,T_ob_eta);
MasValueI3322(4)=max(eig(I3322Op));
%-++ 5
T_ob_eta(:,:,1)=-T_ob_eta(:,:,1);
T_ob_eta(:,:,2)=-T_ob_eta(:,:,2);
T_ob_eta(:,:,3)=-T_ob_eta(:,:,3);
I3322Op=MakeBellOperatorBipartite(I3322_FC,T_ob_eta,T_ob_eta);
MasValueI3322(5)=max(eig(I3322Op));
%-+-
T_ob_eta(:,:,1)=-T_ob_eta(:,:,1);
% T_ob_eta(:,:,2)=-T_ob_eta(:,:,2);
T_ob_eta(:,:,3)=-T_ob_eta(:,:,3);
I3322Op=MakeBellOperatorBipartite(I3322_FC,T_ob_eta,T_ob_eta);
MasValueI3322(6)=max(eig(I3322Op));
%--+
% T_ob_eta(:,:,1)=-T_ob_eta(:,:,1);
T_ob_eta(:,:,2)=-T_ob_eta(:,:,2);
T_ob_eta(:,:,3)=-T_ob_eta(:,:,3);
I3322Op=MakeBellOperatorBipartite(I3322_FC,T_ob_eta,T_ob_eta);
MasValueI3322(7)=max(eig(I3322Op));
%---
% T_ob_eta(:,:,1)=-T_ob_eta(:,:,1);
T_ob_eta(:,:,2)=-T_ob_eta(:,:,2);
T_ob_eta(:,:,3)=-T_ob_eta(:,:,3);
I3322Op=MakeBellOperatorBipartite(I3322_FC,T_ob_eta,T_ob_eta);
MasValueI3322(8)=max(eig(I3322Op));
MasValueI3322=MasValueI3322
if max(MasValueI3322)<=4
    display(['For eta equals to ', num2str(eta)]);
    display(['These measurements do not violate I3322 inequalities (The local bound is 4)']);
else
    display(['For eta equals to ', num2str(eta)]);
    display(['These measurements violate I3322 inequalities (The local bound is 4)']);
end

%CHSH 12
T_ob2_eta=T_ob_eta(:,:,[1 2]);
CHSH=[0 0 0;0 1 1;0 1 -1];
CHSHOp=MakeBellOperatorBipartite(CHSH,T_ob2_eta,T_ob2_eta);
MaxValueCHSH(1)=max(eig(CHSHOp));

%CHSH 13
T_ob2_eta=T_ob_eta(:,:,[1 3]);
CHSH=[0 0 0;0 1 1;0 1 -1];
CHSHOp=MakeBellOperatorBipartite(CHSH,T_ob2_eta,T_ob2_eta);
MaxValueCHSH(2)=max(eig(CHSHOp));

%CHSH 23
T_ob2_eta=T_ob_eta(:,:,[1 3]);
CHSH=[0 0 0;0 1 1;0 1 -1];
CHSHOp=MakeBellOperatorBipartite(CHSH,T_ob2_eta,T_ob2_eta);
MaxValueCHSH(3)=max(eig(CHSHOp));
MaxValueCHSH = MaxValueCHSH

if max(MaxValueCHSH)<=2
    display(['For eta equals to ', num2str(eta)]);
    display(['These measurements do not violate CHSH inequalities (The local bound is 2)']);
else
    display(['For eta equals to ', num2str(eta)]);
    display(['These measurements violate CHSH inequalities (The local bound is 2)']);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% CHANGE VISIBILITY eta %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display(' '); 
display(['           We have analysed the case of visibility ', num2str(eta)]);

eta=.8945;

display(['           We now consider the critial I3322 visibility, which is ', num2str(eta)]);
% display(' ');  

a0=[cos(2*0*pi/3) 0 sin(2*0*pi/3)];
a1=[cos(2*1*pi/3) 0 sin(2*1*pi/3)];
a2=[cos(2*2*pi/3) 0 sin(2*2*pi/3)];

T(:,:,1,1)=Bloch2Rho(a0);
T(:,:,1,2)=Bloch2Rho(a1);
T(:,:,1,3)=Bloch2Rho(a2);
for x=1:3
    T(:,:,2,x)=Id-T(:,:,1,x);
    Teta(:,:,1,x)=eta*T(:,:,1,x)+(1-eta)*Id/2;
    Teta(:,:,2,x)=Id-Teta(:,:,1,x);
end
for x=1:3
    T_ob(:,:,x)=T(:,:,1,x)-T(:,:,2,x);
    T_ob_eta(:,:,x)=eta*T_ob(:,:,x);
end

I3322_FC=[0 1 1 0;-1 1 1 1;-1 1 1 -1;0 1 -1 0];
% I3322_FC(1,1)=-4;


%+++  1
I3322Op=MakeBellOperatorBipartite(I3322_FC,T_ob_eta,T_ob_eta);
MasValueI3322(1)=max(eig(I3322Op));
%++-  2
T_ob_eta(:,:,3)=-T_ob_eta(:,:,3);
I3322Op=MakeBellOperatorBipartite(I3322_FC,T_ob_eta,T_ob_eta);
MasValueI3322(2)=max(eig(I3322Op));
%+-+  3
T_ob_eta(:,:,2)=-T_ob_eta(:,:,2);
T_ob_eta(:,:,3)=-T_ob_eta(:,:,3);
I3322Op=MakeBellOperatorBipartite(I3322_FC,T_ob_eta,T_ob_eta);
MasValueI3322(3)=max(eig(I3322Op));
%+-- 4
T_ob_eta(:,:,2)=-T_ob_eta(:,:,2);
T_ob_eta(:,:,3)=-T_ob_eta(:,:,3);
I3322Op=MakeBellOperatorBipartite(I3322_FC,T_ob_eta,T_ob_eta);
MasValueI3322(4)=max(eig(I3322Op));
%-++ 5
T_ob_eta(:,:,1)=-T_ob_eta(:,:,1);
T_ob_eta(:,:,2)=-T_ob_eta(:,:,2);
T_ob_eta(:,:,3)=-T_ob_eta(:,:,3);
I3322Op=MakeBellOperatorBipartite(I3322_FC,T_ob_eta,T_ob_eta);
MasValueI3322(5)=max(eig(I3322Op));
%-+-
T_ob_eta(:,:,1)=-T_ob_eta(:,:,1);
% T_ob_eta(:,:,2)=-T_ob_eta(:,:,2);
T_ob_eta(:,:,3)=-T_ob_eta(:,:,3);
I3322Op=MakeBellOperatorBipartite(I3322_FC,T_ob_eta,T_ob_eta);
MasValueI3322(6)=max(eig(I3322Op));
%--+
% T_ob_eta(:,:,1)=-T_ob_eta(:,:,1);
T_ob_eta(:,:,2)=-T_ob_eta(:,:,2);
T_ob_eta(:,:,3)=-T_ob_eta(:,:,3);
I3322Op=MakeBellOperatorBipartite(I3322_FC,T_ob_eta,T_ob_eta);
MasValueI3322(7)=max(eig(I3322Op));
%---
% T_ob_eta(:,:,1)=-T_ob_eta(:,:,1);
T_ob_eta(:,:,2)=-T_ob_eta(:,:,2);
T_ob_eta(:,:,3)=-T_ob_eta(:,:,3);
I3322Op=MakeBellOperatorBipartite(I3322_FC,T_ob_eta,T_ob_eta);
MasValueI3322(8)=max(eig(I3322Op));
MasValueI3322=MasValueI3322
if max(MasValueI3322)<=4
    display(['For eta equals to ', num2str(eta)]);
    display(['These measurements do not violate I3322 inequalities (The local bound is 4)']);
else
    display(['For eta equals to ', num2str(eta)]);
    display(['These measurements violate I3322 inequalities (The local bound is 4)']);
end

%CHSH 12
T_ob2_eta=T_ob_eta(:,:,[1 2]);
CHSH=[0 0 0;0 1 1;0 1 -1];
CHSHOp=MakeBellOperatorBipartite(CHSH,T_ob2_eta,T_ob2_eta);
MaxValueCHSH(1)=max(eig(CHSHOp));

%CHSH 13
T_ob2_eta=T_ob_eta(:,:,[1 3]);
CHSH=[0 0 0;0 1 1;0 1 -1];
CHSHOp=MakeBellOperatorBipartite(CHSH,T_ob2_eta,T_ob2_eta);
MaxValueCHSH(2)=max(eig(CHSHOp));

%CHSH 23
T_ob2_eta=T_ob_eta(:,:,[1 3]);
CHSH=[0 0 0;0 1 1;0 1 -1];
CHSHOp=MakeBellOperatorBipartite(CHSH,T_ob2_eta,T_ob2_eta);
MaxValueCHSH(3)=max(eig(CHSHOp));
MaxValueCHSH = MaxValueCHSH

if max(MaxValueCHSH)<=2
    display(['For eta equals to ', num2str(eta)]);
    display(['These measurements do not violate CHSH inequalities (The local bound is 2)']);
else
    display(['For eta equals to ', num2str(eta)]);
    display(['These measurements violate CHSH inequalities (The local bound is 2)']);
end

function [Rho] = Bloch2Rho(BlochVector)
%Function that transforms a Bloch vector into a density matrix
%Input:  BlochVector, the 3-dimensional real vector BlochVector
%Output: Rho, corresponding matrix Rho, if BlochVector nas norm 1, Rho is a density matrix

Rho =1/2*[1+BlochVector(3), BlochVector(1)-sqrt(-1)*BlochVector(2); BlochVector(1)+sqrt(-1)*BlochVector(2), 1-BlochVector(3)];

end

function BellOp=MakeBellOperatorBipartite(F,A,B)
%1 stands for identity
%2 stands for first input
%3 stands for second input

% count the number of inputs
iA=size(A,3);
iB=size(B,3);
d=size(A,1);

%We now shift the operators to set the very first one as identity
% A2(:,:,[2 iA+1])=A
% B2(:,:,[2 iB+1])=B

A2(:,:,1)=eye(d);
B2(:,:,1)=eye(d);
for x=1:iA
    A2(:,:,x+1)=A(:,:,x);
end
for y=1:iB
    B2(:,:,y+1)=B(:,:,y);
end
BellOp=0;
for x=1:iA+1
    for y=1:iB+1
        BellOp=BellOp+F(x,y)*kron(A2(:,:,x),B2(:,:,y));
    end
end
end


