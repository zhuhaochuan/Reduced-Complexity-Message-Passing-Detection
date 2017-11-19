function [J,Z,N0v] = ESTIMATE(N0,N,K,y,Es,H)
N0v = N0/(2*N);
p = sqrt(K*Es);
%Wp = sqrt(N0)*eye(2*N,2*K);
Wp = sqrt(N0v/2) * randn(2*N,2*K);%这里对噪声矩阵的构造 
%Wp = [real(Wp),-imag(Wp);imag(Wp),real(Wp)];
Yp = p*H+Wp;

J = (Yp'*Yp/(N*p^2)-N0v/(p^2)*eye(2*K));
Z = (Yp'*y/(N*p));
