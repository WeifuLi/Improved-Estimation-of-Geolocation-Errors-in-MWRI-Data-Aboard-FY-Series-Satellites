function  [g1]=SurfaceFitting(Data_cubic,Lat1,Data1,Step,Step1)
PadSize=20;
Data_cubic1= padarray(Data_cubic,[PadSize PadSize],'replicate');
%%  model
C=1:1/Step:size(Lat1,1);
Row=1:1/Step1:size(Lat1,2);
[IY,IX] = meshgrid(Row,C);
%%  define A,Z, Kx and Ky
A=zeros(size(Lat1));Z=zeros(size(Lat1));
for ki=1:size(Data1,1)
    for kj=1:size(Data1,2)
        A(IX(ki,kj),IY(ki,kj))=1;
        Z(IX(ki,kj),IY(ki,kj))=Data1(ki,kj);
    end
end
A=padarray(A,[PadSize PadSize]);%% Matlab 数组填充函数
Z=padarray(Z,[PadSize PadSize]);
K=[1 1 1; 1 -8 1; 1 1 1];%  K是对称的, 产生double型的。若不是对称的，则是complex型
sizeA=size(A);
otfF=psf2otf(K,sizeA);
% solution
lambda=1; mu0=0.1; muMax=10^9;
Denormin2 = abs(otfF).^2;
T=mu0+lambda*Denormin2;
Xt=ifft2(mu0*fft2(Data_cubic1)./T);
Cn=0;
while mu0<muMax
    mu0=mu0*2;
    Cn=Cn+1;
    %Fix Xt, solve Zt
    Zt=(A.*Z+mu0*Xt)./(A.*A+mu0);
    %Fix Zt, solve Xt
    Xt=ifft2(mu0*fft2(Zt)./(mu0+lambda*Denormin2));
end
g1=Xt(PadSize+1:end-PadSize,PadSize+1:end-PadSize);
end