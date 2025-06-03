clearvars
clear all
 %%%%%%%%%%%%%%%%%%%%%
%load S10_45_75_90_delta1_2; 
load Sk_deltak1_deltak2_w1_w2_Tabla1_Yang_2016 

lam=[.8,.7,.36,.48];
w1=omega1;
 w2=omega2;
H2=((deltak1).*conj(w2)-(deltak2).*conj(w1))./(w1.*conj(w2)-w2.*conj(w1));
for ii=1:4
 S1=Sk(:,ii);
 delta1=deltak1(ii);
 delta2=deltak2(ii);
 w1=omega1(ii);
 w2=omega2(ii);
 V = abs(w1).*abs(w2).*sin(theta(ii)); % *** Volumen de la celda periodica *******
% 

H1=(pi)/imag(conj(w1)*w2);
H2=((delta1)*conj(w2)-(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
h=(H1);
%h21=imag(H1);
h12=real(H2);
if ii==4
h22=imag(H2);
else
    h22=0;
end
cm=1;
cf=50;
xas=cf/cm;
X=(1-xas)/(1+xas);
     V2=lam(ii); % volumen de la fibra
    r=sqrt(V*V2/pi);  % radio fibra 
    %J=eye(2)+X*r^2*[(H1+h12),(-h22);(-h22),(H1-h12)]
    J=eye(2)+X*r^2*[(h+h12),(-h22);(-h22),(h-h12)];

 %%% usando la notacion de Yang 2016 para la aproximacion
 CENo1(ii,:)=coeficientes_efectivos2f(cm,V2,X,J) ;
 CENo2(ii,:) =ahm2f_Paralelogramo_Cortas_2023_No1(cm,X,V2,S1,r,J);
 CENo3(ii,:) =ahm2f_Paralelogramo_Cortas_2023_No2(cm,X,V2,S1,r,J);
 CENo4(ii,:) =ahm2f_Paralelogramo_Cortas_2023_No3(cm,X,V2,S1,r,J);
%  CENo4(ii,:) =ahm_imperfect_Paralelogramo_Cortas_2023_No4(cm,xas,K,V2,S1,J,r1);
 CENo6(ii,:) =ahm2f_Paralelogramo_Cortas_2023_No5(cm,X,V2,S1,r,J);
 CENo8(ii,:) =conductividad_efectiva_orden_m(cm,X,V2,S1,r,J,7);
 %CENoG(ii,:) = conductividad_efectiva_orden_m(cm,X,V2,S1,r,J,20); 
CENoG(ii,:) = ahm_imperfect_Paralelogramo(cm,cf,V2,V,S1,H1,H2,10);
%  k22yan=cm*(1-(1-3.43759/pi)*X*V2)/(1+(1+3.43759/pi)*X*V2);  
%  k11yan=cm*(1-(1+3.43759/pi)*X*V2)/(1+(1-3.43759/pi)*X*V2);
%  CEyan(ii,:)=[k11yan,k22yan];
 
end



matriz=[ 'cm= ',num2str(cm,5)];
fibra=[ 'cf= ',num2str(cf,5)];
rho=[ '$\rho$','=',num2str(cf/cm,5)];
angulo=['$theta=$ ' num2str(theta*180/pi,3),'$^o$'];
Modulo=[ '|w2|= ',num2str(abs(omega2),5)];
volfibra=[ '$V_2= $',num2str(lam,3)];

disp(Modulo)
disp(angulo)
disp(matriz)
disp(fibra)
disp(volfibra)

textoY=['$\hat{\kappa}_{11}/\kappa_1$';'$\hat{\kappa}_{22}/\kappa_1$';'$\hat{\kappa}_{12}/\kappa_1$'];

Orden=[1:1:4, 6, 8, 10];

Tabla(:,:)=[CENo1(1,1), CENo1(2,1),CENo1(3,1),CENo1(3,2),CENo1(4,1),CENo1(4,2),CENo1(4,3);
            CENo2(1,1), CENo2(2,1),CENo2(3,1),CENo2(3,2),CENo2(4,1),CENo2(4,2),CENo2(4,3);
            CENo3(1,1), CENo3(2,1),CENo3(3,1),CENo3(3,2),CENo3(4,1),CENo3(4,2),CENo3(4,3);
            CENo4(1,1), CENo4(2,1),CENo4(3,1),CENo4(3,2),CENo4(4,1),CENo4(4,2),CENo4(4,3);
            CENo6(1,1), CENo6(2,1),CENo6(3,1),CENo6(3,2),CENo6(4,1),CENo6(4,2),CENo6(4,3);
            CENo8(1,1), CENo8(2,1),CENo8(3,1),CENo8(3,2),CENo8(4,1),CENo8(4,2),CENo8(4,3);
            CENoG(1,1), CENoG(2,1),CENoG(3,1),CENoG(3,2),CENoG(4,1),CENoG(4,2),CENoG(4,3)];  
%%% table 1 Yan 2016
Yan=[7.64407 5.10778 1.66993 3.50961 2.69362 3.79330 0.943634;
     7.64407 6.18007 1.68932 4.02376 2.90209 4.32367 1.22573;
     8.21305 6.29151 1.68987 4.15497 2.94546 4.45468 1.29986;
     8.25787 6.32727 1.68988 4.18703 2.95635 4.48683 1.31857;
     8.25927 6.33552 1.68988 4.19710 2.95970 4.49690 1.32437;
     8.25999 6.33590 1.68988 4.19779 2.95993 4.49759 1.32476;
     8.26001 6.33592 1.68988 4.19784 2.95995 4.49764 1.32479];
 
Table3=[Yan(:,3),Tabla(:,3),Yan(:,4),Tabla(:,4),Yan(:,5),Tabla(:,5),Yan(:,6),Tabla(:,6),Yan(:,7),Tabla(:,7) ] ;      

 disp('O_k       Hex           Sqr       Rect:k11    Rect:k22     Romb:k11    Romb:k22    Romb:k12    ')    
  disp(' -------------------------------------------')
 disp('Table 3')
 disp(num2str([Orden', Table3(:,:)],6))
 
 V2me=[0.90,   0.785,   0.392  0.523]; % percolaci[on reducido
lam=V2me;
 for ii=1:4
 S1=Sk(:,ii);
 delta1=deltak1(ii);
 delta2=deltak2(ii);
 w1=omega1(ii);
 w2=omega2(ii);
 V = abs(w1).*abs(w2).*sin(theta(ii)); % *** Volumen de la celda periodica *******
% 
%H1=(conj(delta1)*conj(w2)-conj(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
H1=(pi)/imag(conj(w1)*w2);
H2=((delta1)*conj(w2)-(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
cm=1;
cf=20;
xas=cf/cm;
X=(1-xas)/(1+xas);
     V2=lam(ii); % volumen de la fibra
    r=sqrt(V*V2/pi);  % radio fibra 
    J=eye(2)+X*r^2*[(h11+h12),(-h22);(-h22),(h11-h12)];
CENoG1(ii,:) = conductividad_efectiva_orden_m(cm,X,V2,S1,r,J,20);
CENoG2(ii,:) = conductividad_efectiva_orden_m(cm,X,V2,S1,r,J,30);
CENoG3(ii,:) = conductividad_efectiva_orden_m(cm,X,V2,S1,r,J,42); 
end

 Vol_per=[ '$V_2= $',num2str(V2me,3)];
 Tabla1(:,:)=[20, CENoG1(1,1),CENoG1(2,1),CENoG1(3,1),CENoG1(3,2),CENoG1(4,1),CENoG1(4,2),CENoG1(4,3);
             30, CENoG2(1,1),CENoG2(2,1),CENoG2(3,1),CENoG2(3,2),CENoG2(4,1),CENoG2(4,2),CENoG2(4,3);
             42, CENoG3(1,1),CENoG3(2,1),CENoG3(3,1),CENoG3(3,2),CENoG3(4,1),CENoG3(4,2),CENoG3(4,3)];

disp(' -------------------------------------------')
 disp(Modulo)
disp(angulo)
disp(matriz)
disp(fibra)
disp(Vol_per)

 disp(' -------------------------------------------')
 disp('Table 4')
 disp(num2str([Tabla1(:,:)],6))                


%%%FUNCIONES
%%% FUNCIONES AUXILIARES


% Funciones auxiliares
%%% orden del sistema No=1
function RR1 =ahm2f_Paralelogramo_Cortas_2023_No1(cm,X,V2,S1,r1,J)

S4=S1(4);S6=S1(6);
L4=[real(S4),-imag(S4);-imag(S4), -real(S4)];
L6=[real(S6),-imag(S6);-imag(S6), -real(S6)];
N1=sqrt(3)*r1^4*L4;
N2=N1;
W=10*r1^6*L6;
AP2=eye(2)+X*W;
PP=-X^2*N1*inv(AP2)*N2;  
JK=J+PP; % Matriz Z
Z=JK;

RR1 =coeficientes_efectivos2f(cm,V2,X,Z);
      
end

%%% orden del sistema No=2
function RR2 =ahm2f_Paralelogramo_Cortas_2023_No2(cm,X,V2,S1,r1,J)
                                                  
S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);
L4=[real(S4),-imag(S4);-imag(S4), -real(S4)];
L6=[real(S6),-imag(S6);-imag(S6), -real(S6)];
L8=[real(S8),-imag(S8);-imag(S8), -real(S8)];
L10=[real(S10),-imag(S10);-imag(S10), -real(S10)];

N11=sqrt(3)*r1^4*L4; N12=sqrt(5)*r1^6*L6;

N1=[N11, N12];
N2=N1';
 
W11=10*r1^6*L6;             W12=7*sqrt(15)*r1^8*L8;              
W21=7*sqrt(15)*r1^8*L8;     W22=126*r1^10*L10;                
 
W=[W11, W12;  W21, W22];

AP2=eye(4)+X*W;
PP=-X^2*N1*inv(AP2)*N2;  
JK=J+PP; % Matriz Z
Z=JK;

RR2 =coeficientes_efectivos2f(cm,V2,X,Z);

end

%%% orden del sistema No=3
function RR3 =ahm2f_Paralelogramo_Cortas_2023_No3(cm,X,V2,S1,r1,J)
              
S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);

L4=[real(S4),-imag(S4);-imag(S4), -real(S4)];
L6=[real(S6),-imag(S6);-imag(S6), -real(S6)];
L8=[real(S8),-imag(S8);-imag(S8), -real(S8)];
L10=[real(S10),-imag(S10);-imag(S10), -real(S10)];
L12=[real(S12),-imag(S12);-imag(S12), -real(S12)];
L14=[real(S14),-imag(S14);-imag(S14), -real(S14)];

N11=sqrt(3)*r1^4*L4; N12=sqrt(5)*r1^6*L6; N13=sqrt(7)*r1^8*L8;
N1=[N11, N12, N13];

N2=N1';
W11=10*r1^6*L6;             W12=7*sqrt(15)*r1^8*L8;      W13=12*sqrt(21)*r1^10*L10;         
W21=7*sqrt(15)*r1^8*L8;     W22=126*r1^10*L10;           W23=66*sqrt(35)*r1^12*L12;        
W31=12*sqrt(21)*r1^10*L10;  W32=66*sqrt(35)*r1^12*L12;   W33=1716*r1^14*L14;                        
W=[W11, W12, W13; W21, W22, W23; W31, W32, W33];

AP2=eye(6)+X*W;
PP=-X^2*N1*inv(AP2)*N2;  
JK=J+PP; % Matriz Z
Z=JK;

RR3 =coeficientes_efectivos2f(cm,V2,X,Z);

end

%%% orden del sistema No=4
function RR3 =ahm2f_Paralelogramo_Cortas_2023_No4(cm,X,V2,S1,r1,J)

S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);S16=S1(16);S18=S1(18);

L4=[real(S4),-imag(S4);-imag(S4), -real(S4)];
L6=[real(S6),-imag(S6);-imag(S6), -real(S6)];
L8=[real(S8),-imag(S8);-imag(S8), -real(S8)];
L10=[real(S10),-imag(S10);-imag(S10), -real(S10)];
L12=[real(S12),-imag(S12);-imag(S12), -real(S12)];
L14=[real(S14),-imag(S14);-imag(S14), -real(S14)];
L16=[real(S16),-imag(S16);-imag(S16), -real(S16)];
L18=[real(S18),-imag(S18);-imag(S18), -real(S18)];

N11=sqrt(3)*r1^4*L4; N12=sqrt(5)*r1^6*L6; N13=sqrt(7)*r1^8*L8; N14=3*r1^10*L10;
N1=[N11, N12, N13, N14];
N2=N1';
W11=10*r1^6*L6;             W12=7*sqrt(15)*r1^8*L8;      W13=12*sqrt(21)*r1^10*L10;    W14=55*sqrt(3)*r1^12*L12;     
W21=7*sqrt(15)*r1^8*L8;     W22=126*r1^10*L10;           W23=66*sqrt(35)*r1^12*L12;    W24=429*sqrt(5)*r1^14*L14;     
W31=12*sqrt(21)*r1^10*L10;  W32=66*sqrt(35)*r1^12*L12;   W33=1716*r1^14*L14;           W34=715*sqrt(63)*r1^16*L16;     
W41=55*sqrt(3)*r1^12*L12;  W42=429*sqrt(5)*r1^14*L14;   W43=715*sqrt(63)*r1^16*L16;   W44=24310*r1^18*L18;  

W=[W11, W12, W13, W14; W21, W22, W23, W24;W31, W32, W33, W34; W41, W42, W43, W44];
    
AP2=eye(8)+X*W;
PP=-X^2*N1*inv(AP2)*N2;  
JK=J+PP; % Matriz Z
Z=JK;

RR3 =coeficientes_efectivos2f(cm,V2,X,Z);

end

%%% orden del sistema No=5
function RR3 =ahm2f_Paralelogramo_Cortas_2023_No5(cm,X,V2,S1,r1,J)
              

S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);S16=S1(16);S18=S1(18);S20=S1(20);S22=S1(22);
L4=[real(S4),-imag(S4);-imag(S4), -real(S4)];
L6=[real(S6),-imag(S6);-imag(S6), -real(S6)];
L8=[real(S8),-imag(S8);-imag(S8), -real(S8)];
L10=[real(S10),-imag(S10);-imag(S10), -real(S10)];
L12=[real(S12),-imag(S12);-imag(S12), -real(S12)];
L14=[real(S14),-imag(S14);-imag(S14), -real(S14)];
L16=[real(S16),-imag(S16);-imag(S16), -real(S16)];
L18=[real(S18),-imag(S18);-imag(S18), -real(S18)];
L20=[real(S20),-imag(S20);-imag(S20), -real(S20)];
L22=[real(S22),-imag(S22);-imag(S22), -real(S22)];

N11=sqrt(3)*r1^4*L4; N12=sqrt(5)*r1^6*L6; N13=sqrt(7)*r1^8*L8; N14=3*r1^10*L10; N15=sqrt(11)*r1^12*L12;

N1=[N11, N12, N13,N14, N15];

N2=N1';

W11=10*r1^6*L6;             W12=7*sqrt(15)*r1^8*L8;      W13=12*sqrt(21)*r1^10*L10;    W14=55*sqrt(3)*r1^12*L12;     W15=26*sqrt(33)*r1^14*L14;
W21=7*sqrt(15)*r1^8*L8;     W22=126*r1^10*L10;           W23=66*sqrt(35)*r1^12*L12;    W24=429*sqrt(5)*r1^14*L14;     W25=273*sqrt(55)*r1^16*L16;
W31=12*sqrt(21)*r1^10*L10;  W32=66*sqrt(35)*r1^12*L12;   W33=1716*r1^14*L14;           W34=2145*sqrt(7)*r1^16*L16;    W35=1768*sqrt(77)*r1^18*L18; 
W41=55*sqrt(3)*r1^12*L12;  W42=429*sqrt(5)*r1^14*L14;   W43=2145*sqrt(7)*r1^16*L16;   W44=24310*r1^18*L18;           W45=25194*sqrt(11)*r1^20*L20; 
W51=26*sqrt(33)*r1^14*L14; W52=273*sqrt(55)*r1^16*L16; W53=1768*sqrt(77)*r1^18*L18; W54=25194*sqrt(11)*r1^20*L20; W55=352716*r1^22*L22;

W=[W11, W12, W13, W14, W15; W21, W22, W23, W24, W25;W31, W32, W33, W34, W35;...
       W41, W42, W43, W44, W45; W51, W52, W53, W54, W55];
   
AP2=eye(10)+X*W;
PP=-X^2*N1*inv(AP2)*N2;  
JK=J+PP; % Matriz Z
Z=JK;

RR3 =coeficientes_efectivos2f(cm,V2,X,Z);

end

%%% orden del sistema No=6
function RR6 =ahm2f_Paralelogramo_Cortas_2023_No6(cm,X,V2,S1,r1,J)

S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);S16=S1(16);S18=S1(18);
S20=S1(20);S22=S1(22);S24=S1(24);S26=S1(26);

L4=[real(S4),-imag(S4);-imag(S4), -real(S4)];
L6=[real(S6),-imag(S6);-imag(S6), -real(S6)];
L8=[real(S8),-imag(S8);-imag(S8), -real(S8)];
L10=[real(S10),-imag(S10);-imag(S10), -real(S10)];
L12=[real(S12),-imag(S12);-imag(S12), -real(S12)];
L14=[real(S14),-imag(S14);-imag(S14), -real(S14)];
L16=[real(S16),-imag(S16);-imag(S16), -real(S16)];
L18=[real(S18),-imag(S18);-imag(S18), -real(S18)];
L20=[real(S20),-imag(S20);-imag(S20), -real(S20)];
L22=[real(S22),-imag(S22);-imag(S22), -real(S22)];
L24=[real(S24),-imag(S24);-imag(S24), -real(S24)];
L26=[real(S26),-imag(S26);-imag(S26), -real(S26)];

N11=sqrt(3)*r1^4*L4; N12=sqrt(5)*r1^6*L6; N13=sqrt(7)*r1^8*L8; N14=3*r1^10*L10; N15=sqrt(11)*r1^12*L12; N16=sqrt(13)*r1^14*L14;

N1=[N11, N12, N13,N14, N15, N16];

N2=N1';

W11=10*r1^6*L6;             W12=7*sqrt(15)*r1^8*L8;      W13=12*sqrt(21)*r1^10*L10;    W14=55*sqrt(3)*r1^12*L12;     W15=26*sqrt(33)*r1^14*L14;
W21=7*sqrt(15)*r1^8*L8;     W22=126*r1^10*L10;           W23=66*sqrt(35)*r1^12*L12;    W24=429*sqrt(5)*r1^14*L14;     W25=273*sqrt(55)*r1^16*L16;
W31=12*sqrt(21)*r1^10*L10;  W32=66*sqrt(35)*r1^12*L12;   W33=1716*r1^14*L14;           W34=2145*sqrt(7)*r1^16*L16;    W35=1768*sqrt(77)*r1^18*L18; 
W41=55*sqrt(3)*r1^12*L12;  W42=429*sqrt(5)*r1^14*L14;   W43=2145*sqrt(7)*r1^16*L16;   W44=24310*r1^18*L18;           W45=25194*sqrt(11)*r1^20*L20; 
W51=26*sqrt(33)*r1^14*L14; W52=273*sqrt(55)*r1^16*L16; W53=1768*sqrt(77)*r1^18*L18; W54=25194*sqrt(11)*r1^20*L20; W55=352716*r1^22*L22;
W61=35*sqrt(39)*r1^16*L16; W62=476*sqrt(65)*r1^18*L18; W63=3876*sqrt(91)*r1^20*L20; W64=22610*sqrt(117)*r1^22*L22; W65=104006*sqrt(143)*r1^24*L24;
W16=35*sqrt(39)*r1^16*L16; W26=476*sqrt(65)*r1^18*L18;  W36=3876*sqrt(91)*r1^20*L20; W46=22610*sqrt(117)*r1^22*L22; W56=104006*sqrt(143)*r1^24*L24;  
W66=5200300*r1^26*L26;

W=[W11, W12, W13, W14, W15, W16; W21, W22, W23, W24, W25, W26; W31, W32, W33, W34, W35, W36;...
   W41, W42, W43, W44, W45, W46; W51, W52, W53, W54, W55, W56; W61, W62, W63, W64, W65, W66];

AP2=eye(12)+X*W;
PP=-X^2*N1*inv(AP2)*N2;  
JK=J+PP; % Matriz Z
Z=JK;

RR6 =coeficientes_efectivos2f(cm,V2,X,Z);

end


function [kappa_eff] = conductividad_efectiva_orden_m(cm,X,V2,Sk,R,J,m)
                    
    % Determine maximum L matrix needed
    max_L_needed = 4*(m) + 2; % Maximum L index needed: 2(m)+2
    
    % Initialize L matrices
    L = cell(1, max_L_needed);
    for k = 1:max_L_needed
        if k <= length(Sk)
            Sk_real = real(Sk(k));
            Sk_imag = imag(Sk(k));
            L{k} = [Sk_real, -Sk_imag; 
                    -Sk_imag, -Sk_real];
        else
            % If Sk not provided for higher orders, set to identity
            warning('Sk(%d) not provided, using identity matrix', k);
            L{k} = eye(2);
        end
    end
    
    % Construct N1 matrix (2 x 2*(m))
    N1_blocks = cell(1, m);
    for j = 1:m
        L_index = 2*j + 2;
        if L_index <= length(L)
            N1_blocks{j} = sqrt(2*j+1) * R^(2*j+2) * L{L_index};
        else
            error('Required L matrix (L_%d) not available', L_index);
        end
    end
    N1 = [N1_blocks{:}];
    
    N2 = (N1');
    
    % Construct W matrix (2*(m) x 2*(m))
    W = zeros(2*(m));
    for s = 1:m
        k = 2*s + 1;
        for t = 1:m
            p = 2*t + 1;
            L_index = k + p;
            
            % Verify we have all required components
            if L_index > length(L)
                error('Required L matrix (L_%d) not available for W construction', L_index);
            end
                        
            coeff = factorial(k + p - 1)*sqrt(k * p) / (factorial(k) * factorial(p));
            W_block = coeff * R^(k + p) * L{k + p} ;
            
            % Place the block in the correct position
            row_range = (2*s-1):(2*s);
            col_range = (2*t-1):(2*t);
            W(row_range, col_range) = W_block;
        end
    end
    
    AP2=eye(2*(m)) + X*W;
    PP=-X^2*N1*inv(AP2)*N2;  
    JK=J+PP; % Matriz Z
    Z=JK;
   
    % Calculate effective conductivity coefficients
    kappa_eff =coeficientes_efectivos2f(cm,V2,X,Z);
 
end

 
 function RR=ahm_imperfect_Paralelogramo(cm,cf,V2,V,S1,H1,H2,No)  
 
r = sqrt(V*V2/pi);  % *** Radio de la fibra *******************************
xas=cf/cm;
S=S1;
%Matriz Wnk
for n=1:2:91
   for k=1:2:91

Wnk(n,k)=(S(n+k))*(r^(n+k))*((factorial(n+k-1)/(factorial(n-1)*factorial(k-1)))/...
               sqrt(n*k));

   end
end
RWnk=real(Wnk);
IWnk=imag(Wnk);
di=eye(length(S)); 

% Para hallar los Bp

for s=1:2:90,

       Bp(s)=((1-xas))/((xas+1));
       
end
% %************************************************************************

%      %*************** MATRIZ del Sistema Problema 13 y 23 . Problema antiplano
  for t=1:No,
   for s=1:No ;

   AP(:,:,t,s)=di(t,s)*[1,0;0,1]+Bp(2*t+1)*[RWnk(2*s+1,2*t+1),-IWnk(2*s+1,2*t+1);...
                                        -IWnk(2*s+1,2*t+1), -RWnk(2*s+1,2*t+1)];
  
   end
  end
  
  for t=1:No
    NP(:,:,t)=Bp(2*t+1)*[RWnk(1,2*t+1),-IWnk(1,2*t+1);...
                 -IWnk(1,2*t+1), -RWnk(1,2*t+1)];
    NPa(:,:,t)=[RWnk(2*t+1,1),-IWnk(2*t+1,1);...
                 -IWnk(2*t+1,1), -RWnk(2*t+1,1)];
    
 end
 
 AP1=zeros(2*No,2,No);
 for s=1:No
 for t=1:No
     AP1(1+2*(t-1):2+2*(t-1),:,s)=AP(:,:,t,s);
 end
 end
 
 AP2=zeros(2*No,2*No);
 for t=1:No
  
  AP2(:,1+2*(t-1):2+2*(t-1))=AP1(:,:,t);
  end

 NP1=zeros(2*No,2);
 NP1a=zeros(2,2*No);
 for t=1:No
   % para armar las matrices
  NP1(1+2*(t-1):2+2*(t-1),:)=NP(:,:,t);
  NP1a(:,1+2*(t-1):2+2*(t-1))=NPa(:,:,t);
  end 
 
  
  PP=-Bp(1)*NP1a*inv(AP2)*NP1; 
  
 

h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);

% J=[1+Bp(1)*r^2*(h11+h12),Bp(1)*r^2*(h21-h22);...
%      -Bp(1)*r^2*(h21+h22),1+Bp(1)*r^2*(h11-h12)]; % cortas
 
 J=eye(2)+Bp(1)*r^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)];

JK=J+PP;
Z=JK;
DtZ=det(Z);

%********************************************************************************************
% FORMULA DE LOS COEFICIENTES EFECTIVOS (PROBLEMAS Antiplano)
%********************************************************************************************

c55=cm*(1-2*V2*Bp(1)*Z(2,2)/DtZ);
  c45=cm*(2*V2*Bp(1)*Z(2,1)/DtZ);
  c54=cm*(2*V2*Bp(1)*Z(1,2)/DtZ);
c44=cm*(1-2*V2*Bp(1)*Z(1,1)/DtZ);
   

RR=[c55,c44,c54,c45];

 end

%**************************************************************************
% FORMULA DE LOS COEFICIENTES EFECTIVOS (PROBLEMAS Antiplano)
%**************************************************************************
function RR =coeficientes_efectivos2f(cm,V2,X,Z) 

DtZ=det(Z);
c55=cm*(1-2*V2*X*Z(2,2)/DtZ);
c45=cm*(2*V2*X*Z(2,1)/DtZ);
c54=cm*(2*V2*X*Z(1,2)/DtZ);
c44=cm*(1-2*V2*X*Z(1,1)/DtZ);
        
RR=[c55,c44,c54,c45];
end

