clc; clear all; close all
mapa=gray(256);
load ('Mesencefalo.mat')

A = volumen (:,:,18);
A = (A - min(A(:)))*255./(max(A(:) - min(A(:))));
A = double(A);
A = filter2(fspecial('average'),A);
se = strel('disk',22);
A = imtophat(A,se);

[m,n] = size (A);
phiX = zeros(m,n);
phiY = zeros(m,n);
phiXY = zeros(m,n);
phiXX = zeros(m,n);
phiYY= zeros(m,n);

C = imagesc(A); colormap(mapa)
E = imellipse(gca,[190.511520737327 261.959064327489 61.2211981566818 53.7076023391792]);
phi = createMask(E, C);
phi = double(phi);

for i=1:m
    for j=1:n
        if phi(i,j)==0;
            phi(i,j)=-1;
        else
            phi(i,j)=1;
        end
    end
end

for x=1:1100
    epsilon=1.0;
    a=0.1;
    lambda1=1;
    lambda2=3;
    mu=0.2;
    v=0;
    t = 10.0;
    dt=0.2;
    Dphi = (epsilon/pi)./(epsilon^2+phi.^2);
    in = find(phi(:) >= 0);
    out = find(phi(:) < 0);
    cIn = mean(A(in));
    cOut = mean(A(out));
    for i=2:m-1
        for j=2:n-1
            phiX(i,j)= (phi(i+1,j) - phi(i-1,j))/2;
            phiY(i,j)= (phi(i,j+1) - phi(i,j-1))/2;
            phiXY(i,j)= (phi(i+1,j+1) - phi(i-1,j+1) - phi(i+1,j-1) + phi(i-1,j-1))/4;
            phiXX(i,j)= phi(i+1,j) + phi(i-1,j) - 2*phi(i,j);
            phiYY(i,j)= phi(i,j+1) + phi(i,j-1) - 2*phi(i,j);
        end
    end
    num = (phiX.^2.*phiYY) + (phiY.^2.*phiXX) - (2.*phiX.*phiY.*phiXY);
    den = ((phiX.^2 + phiY.^2).^(3/2)) + a;
    k = num./den;
    fuerza = -Dphi.*(mu*k-v-(lambda1*(A - cIn).^2) + (lambda2*(A - cOut).^2));
    P=0.0*(4*del2(phi) - k);
    phi = phi + ((0.02)*fuerza + P);
    figure(1)
    imagesc(A); colormap(mapa); hold on;
    contour(phi, [0 0], 'g','LineWidth',1);
    hold off; title([num2str(x) ' Iteraciones']); drawnow;
end
