clearvars

 %%%%%%%%%%%%%%%%%%%%%
%load S10_45_75_90_delta1_2; 
load Sk_deltak1_deltak2_w1_w2_Tabla1_Yang_2016

lam=[.8,.7,.36,.48];
for ii=1:4
 S1=Sk(:,ii);
 delta1=deltak1(ii);
 delta2=deltak2(ii);
 w1=omega1(ii);
 w2=omega2(ii);
 V = abs(w1).*abs(w2).*sin(theta(ii)); % *** Volumen de la celda periodica *******
% 
H1=(conj(delta1)*conj(w2)-conj(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
H2=((delta1)*conj(w2)-(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
cm=1;
cf=50;
xas=cf/cm;
X=(1-xas)/(1+xas);
     V2=lam(ii); % volumen de la fibra
    r=sqrt(V*V2/pi);  % radio fibra 
    J=eye(2)+X*r^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)];

 %%% usando la notacion de Yang 2016 para la aproximacion
 CENo1(ii,:)=coeficientes_efectivos2f(cm,V2,X,J) ;
 CENo2(ii,:) =ahm2f_Paralelogramo_Cortas_2023_No1(cm,X,V2,S1,r,J);
 CENo3(ii,:) =ahm2f_Paralelogramo_Cortas_2023_No2(cm,X,V2,S1,r,J);
 CENo4(ii,:) =ahm2f_Paralelogramo_Cortas_2023_No3(cm,X,V2,S1,r,J);
%  CENo4(ii,:) =ahm_imperfect_Paralelogramo_Cortas_2023_No4(cm,xas,K,V2,S1,J,r1);
 CENo6(ii,:) =ahm2f_Paralelogramo_Cortas_2023_No5(cm,X,V2,S1,r,J);
 CENo8(ii,:) =conductividad_efectiva_orden_m(cm,X,V2,S1,r,J,7);
 CENoG(ii,:) = conductividad_efectiva_orden_m(cm,X,V2,S1,r,J,9); 
                 
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


textoY=['$\hat{\kappa}_{11}/\kappa_1$';'$\hat{\kappa}_{22}/\kappa_1$';'$\hat{\kappa}_{12}/\kappa_1$'];

Orden=[1:1:4, 6, 8, 10];
Zem=[[26], 8.2603 6.3359 1.8840 4.1978 2.9692 4.4972 1.3256] ; %%%% \cite{ZEMLYANOVA2023},
 Yan=[[25],8.2600 6.3359 1.6899 4.1978 2.9600 4.4977 1.3248] ;  %%%% \cite{Yan2016}

Tabla(:,:)=[CENo1(1,1), CENo1(2,1),CENo1(3,1),CENo1(3,2),CENo1(4,1),CENo1(4,2),CENo1(4,3);
            CENo2(1,1), CENo2(2,1),CENo2(3,1),CENo2(3,2),CENo2(4,1),CENo2(4,2),CENo2(4,3);
            CENo3(1,1), CENo3(2,1),CENo3(3,1),CENo3(3,2),CENo3(4,1),CENo3(4,2),CENo3(4,3);
            CENo4(1,1), CENo4(2,1),CENo4(3,1),CENo4(3,2),CENo4(4,1),CENo4(4,2),CENo4(4,3);
            CENo6(1,1), CENo6(2,1),CENo6(3,1),CENo6(3,2),CENo6(4,1),CENo6(4,2),CENo6(4,3);
            CENo8(1,1), CENo8(2,1),CENo8(3,1),CENo8(3,2),CENo8(4,1),CENo8(4,2),CENo8(4,3);
            CENoG(1,1), CENoG(2,1),CENoG(3,1),CENoG(3,2),CENoG(4,1),CENoG(4,2),CENoG(4,3)];         

 disp('O_k       Hex           Sqr       Rect:k11   Rect:k22     Romb:k11    Romb:k22    Romb:k12    ')    
 disp(' -------------------------------------------')
 disp(num2str([Orden', Tabla(:,:); Yan; Zem],6))
 
 
 


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
