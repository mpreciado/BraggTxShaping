function [ R,T,MA,MB ] = simulador_chirpeadas_pruebas( ro )

factor=1;   %Use to prevent temporal aliasing, increasing the temporal window of the simmulation
c=2.9979e8;

N=factor*length(ro);

r=abs(ro);

phi=angle(ro);

t=sqrt(1-(r.^2));

fase=exp(j*pi*([0:N-1]-(N)/2)/(N));

A=[ones(1,N)];
B=[zeros(1,N)];

MA=zeros(length(ro),N);

MB=zeros(length(ro),N);

for(I=1:length(ro))
    phitot=exp(j*phi(I)).*fase;
    A1=A.*phitot+B.*r(I).*phitot;
    B1=A.*r(I).*(phitot.^(-1))+B.*(phitot.^(-1));
    
    A=A1/t(I);
    B=B1/t(I);
    
    MA(I,:)=A;
    
    MB(I,:)=B;
end

R=conj(B)./A;
T=1./A;