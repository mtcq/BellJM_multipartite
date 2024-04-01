%Author: Marco Túlio Quintino, https://github.com/mtcq, mtcq.mm@gmail.com
%Requires: cvx
%Last update: 18/Mar/2024

% This script shows that the critical visibility for a pair of Pauli
% measurements to become JM3 is smaller than eta = 0.8695

% This code starts assuming rho is the GHZ state and eta=0.912;
% Then using the Bell inequality from the witness, we find a better state
% ketMax from it's maximal eigenvalue
%Then, using this new state rho=ketMax*ketMax', we see that the visibility eta2 = 0.8695 leads to a Bell violation

clear all
eta=0.85;
N=3;
d=2;
Id=eye(2);
 a0=[cos(2*0*pi/3) 0 sin(2*0*pi/3)];
 a1=[cos(2*1*pi/3) 0 sin(2*1*pi/3)];
 a2=[cos(2*2*pi/3) 0 sin(2*2*pi/3)];
% a0=[cos(2*0*pi/3) sin(2*0*pi/3) 0];
% a1=[cos(2*1*pi/3) sin(2*1*pi/3) 0];
% a2=[cos(2*2*pi/3) sin(2*2*pi/3) 0];
M(:,:,1,1)=Bloch2Rho(a0);
M(:,:,1,2)=Bloch2Rho(a1);
M(:,:,1,3)=Bloch2Rho(a2);
for x=1:3
    M(:,:,2,x)=Id-M(:,:,1,x);
    Meta(:,:,1,x)=eta*M(:,:,1,x)+(1-eta)*Id/2;
    Meta(:,:,2,x)=Id-Meta(:,:,1,x);
end

ket000=zeros(d^N);
ket000(1)=1;
ket111=zeros(d^N);
ket111(d^N)=1;
ketGHZ=(ket000+sqrt(-1)*ket111)/sqrt(2);
%ketGHZ=(ket000+ket111)/sqrt(2);
rho=ketGHZ*ketGHZ';
% rho=RandomDensityMatrix(d^N,0,1);
for a=1:2
    for b=1:2
        for c=1:2
            for x=1:3
                for y=1:3
                    for z=1:3
                        pQ(a,b,c,x,y,z)=HS_real(rho,kron(kron(Meta(:,:,a,x),Meta(:,:,b,y)),Meta(:,:,c,z)));
                    end
                end
            end
        end
    end
end

[WNR,gammaBell,p_eta,pC,Daxl,Dbyl]=BellLocalTripartiteWNR(pQ);
eta=eta
WNR=WNR
% return
% B=0;
% for a=1:2
%     for b=1:2
%         for c=1:2
%             for x=1:3
%                 for y=1:3
%                     for z=1:3
%                        B=B+gammaBell(a,b,c,x,y,z)*kron(kron(Meta(:,:,a,x),Meta(:,:,b,y)),Meta(:,:,c,z));
%                     end
%                 end
%             end
%         end
%     end
% end

% [V D] = eig(B)
% ketMax=V(:,8);
% ketMax = ketMax/norm(ketMax);
% rho2=ketMax*ketMax';
% 
% eta2=0.8694;
% for x=1:3
%     M(:,:,2,x)=Id-M(:,:,1,x);
%     Meta(:,:,1,x)=eta2*M(:,:,1,x)+(1-eta2)*Id/2;
%     Meta(:,:,2,x)=Id-Meta(:,:,1,x);
% end
% 
% for a=1:2
%     for b=1:2
%         for c=1:2
%             for x=1:3
%                 for y=1:3
%                     for z=1:3
%                         pQ(a,b,c,x,y,z)=HS_real(rho2,kron(kron(Meta(:,:,a,x),Meta(:,:,b,y)),Meta(:,:,c,z)));
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% [WNR2,gammaBell,p_eta,pC,Daxl,Dbyl]=BellLocalTripartiteWNR(pQ);
% eta2=eta2
% WNR2=WNR2


function [Rho] = Bloch2Rho(BlochVector)
%Function that transforms a Bloch vector into a density matrix
%Input:  BlochVector, the 3-dimensional real vector BlochVector
%Output: Rho, corresponding matrix Rho, if BlochVector nas norm 1, Rho is a density matrix

Rho =1/2*[1+BlochVector(3), BlochVector(1)-sqrt(-1)*BlochVector(2); BlochVector(1)+sqrt(-1)*BlochVector(2), 1-BlochVector(3)];

end


function [WNR,gammaBell,p_eta,pC,Daxl,Dbyl]=BellLocalTripartiteWNR(p)

tic;
dS=size(p);
Oa=dS(1);
Ob=dS(2);
Oc=dS(3);
Ia=dS(4);
Ib=dS(5);
Ic=dS(6);

Daxl=Dax_MATRIX(Oa,Ia);
nLA=Oa^Ia;  %Count the deterministic vertices
Dbyl=Dax_MATRIX(Ob,Ib);
nLB=Ob^Ib;  %Count the deterministic vertices

cvx_begin
variable eta
% variable piAB(nLA,nLB)
variable pC(Oc,Ic,nLA,nLB)


expression p_eta
expression aux
dual variable gammaBellCell{Oa,Ob,Oc,Ia,Ib,Ic}

% for c=1:Oc
%     for z=1:Ic
%         for LA=1:nLA
%             for LB=1:nLB
%                 %                 piAB(LA,LB)>=0
%                 pC(c,z,LA,LB)>=0;
%             end
%         end
%     end
% end
pC(:)>=0;

for a=1:Oa
    for b=1:Ob
        for c=1:Oc
            for x=1:Ia
                for y=1:Ib
                    for z=1:Ic
                        p_eta(a,b,c,x,y,z)=eta*p(a,b,c,x,y,z) + (1-eta)/(Oa*Ob*Oc);
                        aux=0;
                        for LA=1:nLA
                            for LB=1:nLB
                                aux=aux+pC(c,z,LA,LB)*Daxl(a,x,LA)*Dbyl(b,y,LB);
                            end
                        end
                        gammaBellCell{a,b,c,x,y,z} : p_eta(a,b,c,x,y,z)==aux ;
                    end
                end
            end
        end
    end
end

maximise eta
eta<=10;

cvx_end

WNR=eta;
for a=1:Oa
    for b=1:Ob
        for c=1:Oc
            for x=1:Ia
                for y=1:Ib
                    for z=1:Ic
                        gammaBell(a,b,c,x,y,z)=gammaBellCell{a,b,c,x,y,z};
                    end
                end
            end
        end
    end
end

toc
end

function Dax = Dax_MATRIX(Oa,Ia)

L=Oa^Ia; %Number of deterministic strategies
Dax=cell(Oa,Ia,L);

nA=Ia;
mA=Oa;
Ndet = mA^nA;
Idm = eye(mA);
SingleParty = zeros(nA*mA,Ndet);
multf = linspace(nA-1,0,nA); %used for 'binary' filling of arrays
multf = mA.^multf;

% odo stands for odometer
odo_n = nA;             % number of for loops
odo_k = mA*ones(1,nA);  % number of iterations inside each for loop in order
odo_inc = ones(1,nA);   % increment value for each for loop
odo_sp = zeros(1,nA);   %starting position of odometer

odo = odo_sp;
exitflat = 0;
while exitflat < 1
    
    for i = odo_n:-1:1
        if odo(i)+odo_inc(i) <= odo_sp+(odo_k(i)-1)*odo_inc(i)
            odo(i) = odo(i) + odo_inc(i);
            break
        else
            odo(i) = odo_sp(i);
        end
        
        if i == 1
            exitflat = 1;
        end
    end
    
    
    temp = [];
    for j = 1:nA
        temp = [temp; Idm(:,odo(j)+1)];
    end
    SingleParty(:,1+multf*odo') = temp;
    
end

Dax=reshape(SingleParty,Oa,Ia,L);

end

function out=HS_real(A,B)
%Function that evaluates the real Hilbert-Schmidt inner product between two matrices A,B  out=real(real(trace(A'*B))
%This is useful to reduce complexity and reduce numerical imprecision
%A=A', B=B', HS(A,B) is always real, hence this function prevents erros

%Input: A pair or matrices A and B
%Output: real(trace(A'*B))

out=real(A(:)'*B(:));
end