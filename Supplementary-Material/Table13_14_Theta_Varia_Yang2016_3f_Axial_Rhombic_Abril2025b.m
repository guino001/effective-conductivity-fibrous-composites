clearvars

load Sk_deltak1_deltak2_w1_w2_Tabla3_Yang_2016
r=[1,1,sqrt(5)/2, 2*sqrt(2)/(1+sqrt(3))];
theta=[pi/3,pi/2,atan(2),5*pi/12];
omega1=[1,1,1,1];
omega2=r.*exp(1i*theta);
%%%%Percolacion
V1 = abs(omega1).*abs(omega2).*sin(theta); % *** Volumen de la celda periodica *******
d=sqrt((abs(omega1)).^2+(abs(omega2)).^2-2.*abs(omega1).*abs(omega2).*cos(theta));
rm=1/2.*min([abs(omega1);abs(omega2);d]',[],2);
rm=rm';
V2m=pi.*rm.^2./V1; % percolacion

V2masV3=[.8,.7,.7,.7];
xx=0.1; % t/R2
cm=1;
cmf=990.5;
cf=1/100;
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
h12=real(H2);
h22=imag(H2);
x1=cmf/cm;
x2=cf/cm;
% *************************************************************************
V3=V2masV3(ii)/(1+xx).^2; % volumen de fibra interna
V2=V2masV3(ii)-V3; % volumen del anillo
    r1=sqrt(V*V2masV3(ii)/pi);  % radio fibra mas mesofase
    r2=sqrt(V*V3/pi);  % radio fibra mas fibra interior
    t=xx*r2; %% espesor de la mesofase
X1=((1-x1)*(x1+x2)*(V2+V3)+(1+x1)*(x1-x2)*(V3))/...
             ((1+x1)*(x1+x2)*(V2+V3)+(1-x1)*(x1-x2)*(V3));
 J=eye(2)+X1*r1^2*[(h+h12),(-h22);(-h22),(h-h12)];
 
 CENo0(ii,:)=coeficientes_efectivos(x1,x2,V3,V2,J); 
 CENo1(ii,:) =ahm3f_Paralelogramo_Cortas_2025_No1(x1,x2,V3,V2,S1,r1,J);
 CENo2(ii,:) =ahm3f_Paralelogramo_Cortas_2025_No2(x1,x2,V3,V2,S1,r1,J); 
 CENo3(ii,:) =ahm3f_Paralelogramo_Cortas_2025_No3(x1,x2,V3,V2,S1,r1,J);
 CENo4(ii,:) =ahm3f_Paralelogramo_Cortas_2025_No5(x1,x2,V3,V2,S1,r1,J);
 CENo5(ii,:) = conductividad_efectiva_orden_m(cm,cf,cmf,V3,V2,V,S1,r1,J,7);
 CENoG(ii,:) = conductividad_efectiva_orden_m(cm,cf,cmf,V3,V2,V,S1,r1,J,9);
                 
end


V2masV3a=[ 0.906899,0.785398,0.785398,0.785398];
xx=0.1; % t/R2
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
h12=real(H2);
h22=imag(H2);

x1=cmf/cm;
x2=cf/cm;
% *************************************************************************
V3=V2masV3a(ii)/(1+xx).^2; % volumen de fibra interna
V2=V2masV3a(ii)-V3; % volumen del anillo
    r1=sqrt(V*V2masV3a(ii)/pi);  % radio fibra mas mesofase
    r2=sqrt(V*V3/pi);  % radio fibra mas fibra interior
    t=xx*r2; %% espesor de la mesofase
X1=((1-x1)*(x1+x2)*(V2+V3)+(1+x1)*(x1-x2)*(V3))/...
             ((1+x1)*(x1+x2)*(V2+V3)+(1-x1)*(x1-x2)*(V3));
 J=eye(2)+X1*r1^2*[(h+h12),(-h22);(-h22),(h-h12)];
 
 CENoG1(ii,:) = conductividad_efectiva_orden_m(cm,cf,cmf,V3,V2,V,S1,r1,J,22);
 CENoG2(ii,:) = conductividad_efectiva_orden_m(cm,cf,cmf,V3,V2,V,S1,r1,J,32); 
 CENoG3(ii,:) = conductividad_efectiva_orden_m(cm,cf,cmf,V3,V2,V,S1,r1,J,42);  
end

matriz=[ 'cm= ',num2str(cm,5)];
fibra=[ 'cf= ',num2str(cf,5)];
mesofase=[ 'cmf= ',num2str(cmf,5)];
rho1=[ '$\rho_1$','=',num2str(cmf/cm,5),];
rho2=[ '$\rho_2$','=',num2str(cf/cm,5)];
angulo=['$theta=$ ' num2str(theta*180/pi,3),'$^o$'];
Modulo=[ '|w2|= ',num2str(abs(omega2),5)];
volmesofase=[ '$V_2= $',num2str(V2,3)];
disp(Modulo)
disp(angulo)
disp(matriz)
disp(fibra)
disp(mesofase)

textoY=['$\hat{\kappa}_{11}/\kappa_1$';'$\hat{\kappa}_{22}/\kappa_1$';'$\hat{\kappa}_{12}/\kappa_1$'];


Orden=[1:1:4,6,8, 10];
Zemi=[ 9.0532 6.8758 7.3080 4.7771 7.1407 5.5385 0.4449];
Yan =[8.99259 6.83800 7.26517 4.76177 7.09965 5.51602 0.440128];

Yan=[ 8.22428 5.35492 6.47373 4.61585 5.97169 4.94750 0.452284
      8.22428 6.63206 6.90262 4.68798 6.80287 5.44000 0.342860
      8.92889 6.77930 7.22368 4.76057 7.05683 5.50599 0.454821
      8.98960 6.82609 7.25461 4.76175 7.08825 5.51496 0.438090
      8.99155 6.83743 7.26463 4.76176 7.09915 5.51601 0.440132
      8.99255 6.83797 7.26515 4.76177 7.09963 5.51602 0.440131
      8.99259 6.83800 7.26517 4.76177 7.09965 5.51602 0.440128];

Tabla(:,:)=[CENo0(1,1), CENo0(2,1),CENo0(3,1),CENo0(3,2),CENo0(4,1),CENo0(4,2),CENo0(4,3);
            CENo1(1,1), CENo1(2,1),CENo1(3,1),CENo1(3,2),CENo1(4,1),CENo1(4,2),CENo1(4,3);
            CENo2(1,1), CENo2(2,1),CENo2(3,1),CENo2(3,2),CENo2(4,1),CENo2(4,2),CENo2(4,3);
            CENo3(1,1), CENo3(2,1),CENo3(3,1),CENo3(3,2),CENo3(4,1),CENo3(4,2),CENo3(4,3);
            CENo4(1,1), CENo4(2,1),CENo4(3,1),CENo4(3,2),CENo4(4,1),CENo4(4,2),CENo4(4,3);
            CENo5(1,1), CENo5(2,1),CENo5(3,1),CENo5(3,2),CENo5(4,1),CENo5(4,2),CENo5(4,3);
            CENoG(1,1), CENoG(2,1),CENoG(3,1),CENoG(3,2),CENoG(4,1),CENoG(4,2),CENoG(4,3)];    
   Table13=[Yan(:,3),Tabla(:,3),Yan(:,4),Tabla(:,4),Yan(:,5),Tabla(:,5),Yan(:,6),Tabla(:,6),Yan(:,7),Tabla(:,7) ] ;

 disp('O_k       Hex           Sqr       Medio:k11   Medio:k22     Romb:k11    Romb:k22    Romb:k12    ')    
 disp(' -------------------------------------------')
 disp('Table 13')
 disp(num2str([Orden', Table13(:,:)],6))
 Vol_per=[ '$V_f= $',num2str(V2masV3a,7)];
 
 
 Tabla1(:,:)=[22, CENoG1(1,1),CENoG1(2,1),CENoG1(3,1),CENoG1(3,2),CENoG1(4,1),CENoG1(4,2),CENoG1(4,3);
             32, CENoG2(1,1),CENoG2(2,1),CENoG2(3,1),CENoG2(3,2),CENoG2(4,1),CENoG2(4,2),CENoG2(4,3);
             42, CENoG3(1,1),CENoG3(2,1),CENoG3(3,1),CENoG3(3,2),CENoG3(4,1),CENoG3(4,2),CENoG3(4,3)];
                
 disp(Modulo)
disp(angulo)
disp(matriz)
disp(fibra)
disp(mesofase)
disp(Vol_per)

 disp(' -------------------------------------------')
 disp('Table 14')
 disp(num2str([Tabla1(:,:)],6))
 
 


%%%FUNCIONES

%% orden del sistema No=1
function RR1 =ahm3f_Paralelogramo_Cortas_2025_No1(x1,x2,V3,V2,S1,r1,J)

S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);S16=S1(16);S18=S1(18);
S20=S1(20);S22=S1(22);S24=S1(24);
X1=((1-x1)*(x1+x2)*(V2+V3)+(1+x1)*(x1-x2)*(V3))/...
             ((1+x1)*(x1+x2)*(V2+V3)+(1-x1)*(x1-x2)*(V3));
X3=((1-x1)*(x1+x2)*(V2+V3)^3+(1+x1)*(x1-x2)*(V3)^3)/...
             ((1+x1)*(x1+x2)*(V2+V3)^3+(1-x1)*(x1-x2)*(V3)^3);
X5=((1-x1)*(x1+x2)*(V2+V3)^5+(1+x1)*(x1-x2)*(V3)^5)/...
             ((1+x1)*(x1+x2)*(V2+V3)^5+(1-x1)*(x1-x2)*(V3)^5);
X7=((1-x1)*(x1+x2)*(V2+V3)^7+(1+x1)*(x1-x2)*(V3)^7)/...
             ((1+x1)*(x1+x2)*(V2+V3)^7+(1-x1)*(x1-x2)*(V3)^7);
X9=((1-x1)*(x1+x2)*(V2+V3)^9+(1+x1)*(x1-x2)*(V3)^9)/...
             ((1+x1)*(x1+x2)*(V2+V3)^9+(1-x1)*(x1-x2)*(V3)^9);
 
X11=((1-x1)*(x1+x2)*(V2+V3)^11+(1+x1)*(x1-x2)*(V3)^11)/...
             ((1+x1)*(x1+x2)*(V2+V3)^11+(1-x1)*(x1-x2)*(V3)^11);
        
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

N11=sqrt(3)*r1^4*L4;

N1=[N11];

N2=[X3*N11];

W11=10*X3*r1^6*L6;                             
                     

W=[W11];

AP2=eye(2)+W;
PP=-X1*N1*inv(AP2)*N2;  
JK=J+PP; % Matriz Z
% disp('corto');
Z=JK;

RR1 =coeficientes_efectivos(x1,x2,V3,V2,Z);

end



%%% orden del sistema No=2
function RR2 =ahm3f_Paralelogramo_Cortas_2025_No2(x1,x2,V3,V2,S1,r1,J)

S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);S16=S1(16);S18=S1(18);
S20=S1(20);S22=S1(22);S24=S1(24);
X1=((1-x1)*(x1+x2)*(V2+V3)+(1+x1)*(x1-x2)*(V3))/...
             ((1+x1)*(x1+x2)*(V2+V3)+(1-x1)*(x1-x2)*(V3));
X3=((1-x1)*(x1+x2)*(V2+V3)^3+(1+x1)*(x1-x2)*(V3)^3)/...
             ((1+x1)*(x1+x2)*(V2+V3)^3+(1-x1)*(x1-x2)*(V3)^3);
X5=((1-x1)*(x1+x2)*(V2+V3)^5+(1+x1)*(x1-x2)*(V3)^5)/...
             ((1+x1)*(x1+x2)*(V2+V3)^5+(1-x1)*(x1-x2)*(V3)^5);
X7=((1-x1)*(x1+x2)*(V2+V3)^7+(1+x1)*(x1-x2)*(V3)^7)/...
             ((1+x1)*(x1+x2)*(V2+V3)^7+(1-x1)*(x1-x2)*(V3)^7);
X9=((1-x1)*(x1+x2)*(V2+V3)^9+(1+x1)*(x1-x2)*(V3)^9)/...
             ((1+x1)*(x1+x2)*(V2+V3)^9+(1-x1)*(x1-x2)*(V3)^9);
 
X11=((1-x1)*(x1+x2)*(V2+V3)^11+(1+x1)*(x1-x2)*(V3)^11)/...
             ((1+x1)*(x1+x2)*(V2+V3)^11+(1-x1)*(x1-x2)*(V3)^11);
        
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

N11=sqrt(3)*r1^4*L4; N12=sqrt(5)*r1^6*L6;

N1=[N11, N12];

N2=[X3*N11; X5*N12];

W11=10*X3*r1^6*L6;             W12=7*sqrt(15)*X3*r1^8*L8;             
W21=7*sqrt(15)*X5*r1^8*L8;     W22=126*X5*r1^10*L10;                
                     

W=[W11, W12; W21, W22];

AP2=eye(4)+W;
PP=-X1*N1*inv(AP2)*N2;  
JK=J+PP; % Matriz Z
% disp('corto');
Z=JK;

RR2 =coeficientes_efectivos(x1,x2,V3,V2,Z);

end


%%% orden del sistema No=3
function RR3 =ahm3f_Paralelogramo_Cortas_2025_No3(x1,x2,V3,V2,S1,r1,J)

S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);S16=S1(16);S18=S1(18);
S20=S1(20);S22=S1(22);S24=S1(24);
X1=((1-x1)*(x1+x2)*(V2+V3)+(1+x1)*(x1-x2)*(V3))/...
             ((1+x1)*(x1+x2)*(V2+V3)+(1-x1)*(x1-x2)*(V3));
X3=((1-x1)*(x1+x2)*(V2+V3)^3+(1+x1)*(x1-x2)*(V3)^3)/...
             ((1+x1)*(x1+x2)*(V2+V3)^3+(1-x1)*(x1-x2)*(V3)^3);
X5=((1-x1)*(x1+x2)*(V2+V3)^5+(1+x1)*(x1-x2)*(V3)^5)/...
             ((1+x1)*(x1+x2)*(V2+V3)^5+(1-x1)*(x1-x2)*(V3)^5);
X7=((1-x1)*(x1+x2)*(V2+V3)^7+(1+x1)*(x1-x2)*(V3)^7)/...
             ((1+x1)*(x1+x2)*(V2+V3)^7+(1-x1)*(x1-x2)*(V3)^7);
X9=((1-x1)*(x1+x2)*(V2+V3)^9+(1+x1)*(x1-x2)*(V3)^9)/...
             ((1+x1)*(x1+x2)*(V2+V3)^9+(1-x1)*(x1-x2)*(V3)^9);
 
X11=((1-x1)*(x1+x2)*(V2+V3)^11+(1+x1)*(x1-x2)*(V3)^11)/...
             ((1+x1)*(x1+x2)*(V2+V3)^11+(1-x1)*(x1-x2)*(V3)^11);
        
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

N11=sqrt(3)*r1^4*L4; N12=sqrt(5)*r1^6*L6; N13=sqrt(7)*r1^8*L8;

N1=[N11, N12, N13];

N2=[X3*N11; X5*N12; X7*N13];

W11=10*X3*r1^6*L6;             W12=7*sqrt(15)*X3*r1^8*L8;      W13=12*sqrt(21)*X3*r1^10*L10;         
W21=7*sqrt(15)*X5*r1^8*L8;     W22=126*X5*r1^10*L10;           W23=66*sqrt(35)*X5*r1^12*L12;        
W31=12*sqrt(21)*X7*r1^10*L10;  W32=66*sqrt(35)*X7*r1^12*L12;   W33=1716*X7*r1^14*L14;                        

W=[W11, W12, W13; W21, W22, W23; W31, W32, W33];

AP2=eye(6)+W;
PP=-X1*N1*inv(AP2)*N2;  
JK=J+PP; % Matriz Z
% disp('corto');
Z=JK;

RR3 =coeficientes_efectivos(x1,x2,V3,V2,Z);

end


%%% orden del sistema No=4
function RR4 =ahm3f_Paralelogramo_Cortas_2025_No4(x1,x2,V3,V2,S1,r1,J)

S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);S16=S1(16);S18=S1(18);
S20=S1(20);S22=S1(22);
X1=((1-x1)*(x1+x2)*(V2+V3)+(1+x1)*(x1-x2)*(V3))/...
             ((1+x1)*(x1+x2)*(V2+V3)+(1-x1)*(x1-x2)*(V3));
X3=((1-x1)*(x1+x2)*(V2+V3)^3+(1+x1)*(x1-x2)*(V3)^3)/...
             ((1+x1)*(x1+x2)*(V2+V3)^3+(1-x1)*(x1-x2)*(V3)^3);
X5=((1-x1)*(x1+x2)*(V2+V3)^5+(1+x1)*(x1-x2)*(V3)^5)/...
             ((1+x1)*(x1+x2)*(V2+V3)^5+(1-x1)*(x1-x2)*(V3)^5);
X7=((1-x1)*(x1+x2)*(V2+V3)^7+(1+x1)*(x1-x2)*(V3)^7)/...
             ((1+x1)*(x1+x2)*(V2+V3)^7+(1-x1)*(x1-x2)*(V3)^7);
X9=((1-x1)*(x1+x2)*(V2+V3)^9+(1+x1)*(x1-x2)*(V3)^9)/...
             ((1+x1)*(x1+x2)*(V2+V3)^9+(1-x1)*(x1-x2)*(V3)^9);
 
        
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

N2=[X3*N11; X5*N12; X7*N13; X9*N14];

W11=10*X3*r1^6*L6;             W12=7*sqrt(15)*X3*r1^8*L8;      W13=12*sqrt(21)*X3*r1^10*L10;    W14=55*sqrt(3)*X3*r1^12*L12;     
W21=7*sqrt(15)*X5*r1^8*L8;     W22=126*X5*r1^10*L10;           W23=66*sqrt(35)*X5*r1^12*L12;    W24=429*sqrt(5)*X5*r1^14*L14;     
W31=12*sqrt(21)*X7*r1^10*L10;  W32=66*sqrt(35)*X7*r1^12*L12;   W33=1716*X7*r1^14*L14;           W34=2145*sqrt(7)*X7*r1^16*L16;     
W41=55*sqrt(3)*X9*r1^12*L12;  W42=429*sqrt(5)*X9*r1^14*L14;   W43=2145*sqrt(7)*X9*r1^16*L16;   W44=24310*X9*r1^18*L18;  

W=[W11, W12, W13, W14; W21, W22, W23, W24;W31, W32, W33, W34; W41, W42, W43, W44];

AP2=eye(8)+W;
PP=-X1*N1*inv(AP2)*N2;  
JK=J+PP; % Matriz Z
% disp('corto');
Z=JK;

RR4 =coeficientes_efectivos(x1,x2,V3,V2,Z);

end


%%% orden del sistema No=5
function RR5 =ahm3f_Paralelogramo_Cortas_2025_No5(x1,x2,V3,V2,S1,r1,J)

S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);S16=S1(16);S18=S1(18);
S20=S1(20);S22=S1(22);S24=S1(24);
X1=((1-x1)*(x1+x2)*(V2+V3)+(1+x1)*(x1-x2)*(V3))/...
             ((1+x1)*(x1+x2)*(V2+V3)+(1-x1)*(x1-x2)*(V3));
X3=((1-x1)*(x1+x2)*(V2+V3)^3+(1+x1)*(x1-x2)*(V3)^3)/...
             ((1+x1)*(x1+x2)*(V2+V3)^3+(1-x1)*(x1-x2)*(V3)^3);
X5=((1-x1)*(x1+x2)*(V2+V3)^5+(1+x1)*(x1-x2)*(V3)^5)/...
             ((1+x1)*(x1+x2)*(V2+V3)^5+(1-x1)*(x1-x2)*(V3)^5);
X7=((1-x1)*(x1+x2)*(V2+V3)^7+(1+x1)*(x1-x2)*(V3)^7)/...
             ((1+x1)*(x1+x2)*(V2+V3)^7+(1-x1)*(x1-x2)*(V3)^7);
X9=((1-x1)*(x1+x2)*(V2+V3)^9+(1+x1)*(x1-x2)*(V3)^9)/...
             ((1+x1)*(x1+x2)*(V2+V3)^9+(1-x1)*(x1-x2)*(V3)^9);
 
X11=((1-x1)*(x1+x2)*(V2+V3)^11+(1+x1)*(x1-x2)*(V3)^11)/...
             ((1+x1)*(x1+x2)*(V2+V3)^11+(1-x1)*(x1-x2)*(V3)^11);
        
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

N2=[X3*N11; X5*N12; X7*N13; X9*N14; X11*N15];

W11=10*X3*r1^6*L6;             W12=7*sqrt(15)*X3*r1^8*L8;      W13=12*sqrt(21)*X3*r1^10*L10;    W14=55*sqrt(3)*X3*r1^12*L12;     W15=26*sqrt(33)*X3*r1^14*L14;
W21=7*sqrt(15)*X5*r1^8*L8;     W22=126*X5*r1^10*L10;           W23=66*sqrt(35)*X5*r1^12*L12;    W24=429*sqrt(5)*X5*r1^14*L14;     W25=273*sqrt(55)*X5*r1^16*L16;
W31=12*sqrt(21)*X7*r1^10*L10;  W32=66*sqrt(35)*X7*r1^12*L12;   W33=1716*X7*r1^14*L14;           W34=2145*sqrt(7)*X7*r1^16*L16;    W35=1768*sqrt(77)*X7*r1^18*L18; 
W41=55*sqrt(3)*X9*r1^12*L12;  W42=429*sqrt(5)*X9*r1^14*L14;   W43=2145*sqrt(7)*X9*r1^16*L16;   W44=24310*X9*r1^18*L18;           W45=25194*sqrt(11)*X9*r1^20*L20; 
W51=26*sqrt(33)*X11*r1^14*L14; W52=273*sqrt(55)*X11*r1^16*L16; W53=1768*sqrt(77)*X11*r1^18*L18; W54=25194*sqrt(11)*X11*r1^20*L20; W55=352716*X11*r1^22*L22;



W=[W11, W12, W13, W14, W15; W21, W22, W23, W24, W25;W31, W32, W33, W34, W35;...
       W41, W42, W43, W44, W45; W51, W52, W53, W54, W55];

AP2=eye(10)+W;
PP=-X1*N1*inv(AP2)*N2;  
JK=J+PP; % Matriz Z
% disp('corto');
Z=JK;

RR5 =coeficientes_efectivos(x1,x2,V3,V2,Z);

end

%%% orden del sistema No=6
function RR6 =ahm3f_Paralelogramo_Cortas_2025_No6(x1,x2,V3,V2,S1,r1,J)

S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);S16=S1(16);S18=S1(18);
S20=S1(20);S22=S1(22);S24=S1(24);S26=S1(26);
X1=((1-x1)*(x1+x2)*(V2+V3)+(1+x1)*(x1-x2)*(V3))/...
             ((1+x1)*(x1+x2)*(V2+V3)+(1-x1)*(x1-x2)*(V3));
X3=((1-x1)*(x1+x2)*(V2+V3)^3+(1+x1)*(x1-x2)*(V3)^3)/...
             ((1+x1)*(x1+x2)*(V2+V3)^3+(1-x1)*(x1-x2)*(V3)^3);
X5=((1-x1)*(x1+x2)*(V2+V3)^5+(1+x1)*(x1-x2)*(V3)^5)/...
             ((1+x1)*(x1+x2)*(V2+V3)^5+(1-x1)*(x1-x2)*(V3)^5);
X7=((1-x1)*(x1+x2)*(V2+V3)^7+(1+x1)*(x1-x2)*(V3)^7)/...
             ((1+x1)*(x1+x2)*(V2+V3)^7+(1-x1)*(x1-x2)*(V3)^7);
X9=((1-x1)*(x1+x2)*(V2+V3)^9+(1+x1)*(x1-x2)*(V3)^9)/...
             ((1+x1)*(x1+x2)*(V2+V3)^9+(1-x1)*(x1-x2)*(V3)^9);
 
X11=((1-x1)*(x1+x2)*(V2+V3)^11+(1+x1)*(x1-x2)*(V3)^11)/...
             ((1+x1)*(x1+x2)*(V2+V3)^11+(1-x1)*(x1-x2)*(V3)^11);
X13=((1-x1)*(x1+x2)*(V2+V3)^13+(1+x1)*(x1-x2)*(V3)^13)/...
             ((1+x1)*(x1+x2)*(V2+V3)^13+(1-x1)*(x1-x2)*(V3)^13);
        
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

N2=[X3*N11; X5*N12; X7*N13; X9*N14; X11*N15; X13*N16];

W11=10*X3*r1^6*L6;             W12=7*sqrt(15)*X3*r1^8*L8;      W13=12*sqrt(21)*X3*r1^10*L10;    W14=55*sqrt(3)*X3*r1^12*L12;     W15=26*sqrt(33)*X3*r1^14*L14;
W21=7*sqrt(15)*X5*r1^8*L8;     W22=126*X5*r1^10*L10;           W23=66*sqrt(35)*X5*r1^12*L12;    W24=429*sqrt(5)*X5*r1^14*L14;     W25=273*sqrt(55)*X5*r1^16*L16;
W31=12*sqrt(21)*X7*r1^10*L10;  W32=66*sqrt(35)*X7*r1^12*L12;   W33=1716*X7*r1^14*L14;           W34=2145*sqrt(7)*X7*r1^16*L16;    W35=1768*sqrt(77)*X7*r1^18*L18; 
W41=55*sqrt(3)*X9*r1^12*L12;  W42=429*sqrt(5)*X9*r1^14*L14;   W43=2145*sqrt(7)*X9*r1^16*L16;   W44=24310*X9*r1^18*L18;           W45=25194*sqrt(11)*X9*r1^20*L20; 
W51=26*sqrt(33)*X11*r1^14*L14; W52=273*sqrt(55)*X11*r1^16*L16; W53=1768*sqrt(77)*X11*r1^18*L18; W54=25194*sqrt(11)*X11*r1^20*L20; W55=352716*X11*r1^22*L22;
W61=35*sqrt(39)*X13*r1^16*L16; W62=476*sqrt(65)*X13*r1^18*L18; W63=3876*sqrt(91)*X13*r1^20*L20; W64=22610*sqrt(117)*X13*r1^22*L22; W65=104006*sqrt(143)*X13*r1^24*L24;
W16=35*sqrt(39)*X3*r1^16*L16; W26=476*sqrt(65)*X5*r1^18*L18;  W36=3876*sqrt(91)*X7*r1^20*L20; W46=22610*sqrt(117)*X9*r1^22*L22; W56=104006*sqrt(143)*X11*r1^24*L24;  
W66=5200300*X13*r1^26*L26;


W=[W11, W12, W13, W14, W15, W16; W21, W22, W23, W24, W25, W26; W31, W32, W33, W34, W35, W36;...
   W41, W42, W43, W44, W45, W46; W51, W52, W53, W54, W55, W56; W61, W62, W63, W64, W65, W66];

AP2=eye(12)+W;
PP=-X1*N1*inv(AP2)*N2;  
JK=J+PP; % Matriz Z
% disp('corto');
Z=JK;

RR6 =coeficientes_efectivos(x1,x2,V3,V2,Z);

end



function [kappa_eff] = conductividad_efectiva_orden_m(kappa1, kappa3, kappa2,V3,V2,V,Sk,R1,J1,m);
                    
   
      V2plusV3=V2+V3;
      %R1=sqrt(V*V2plusV3/pi);  % radio fibra mas mesofase   
    % Calculate rho values
    rho1 = kappa2 / kappa1;
    rho2 = kappa3 / kappa1;
    
    % Calculate chi values up to required order
    max_p = 6*(m) + 1; % We need up to chi_{2(m)+1}
    chi = zeros(1, max_p);
    for p = 1:2:max_p
        numerator = (1-rho1)*(rho1+rho2)*(V2+V3)^p + (1+rho1)*(rho1-rho2)*V3^p;
        denominator = (1+rho1)*(rho1+rho2)*(V2+V3)^p + (1-rho1)*(rho1-rho2)*V3^p;
        chi(p) = numerator / denominator;
    end
    
    
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
            N1_blocks{j} = sqrt(2*j+1) * R1^(2*j+2) * L{L_index};
        else
            error('Required L matrix (L_%d) not available', L_index);
        end
    end
    N1 = [N1_blocks{:}];
    
    % Construct N2 matrix (2*(m) x 2)
    N2_blocks = cell(m, 1);
    for j = 1:m
        chi_index = 2*j + 1;
        if chi_index <= length(chi)
            N2_blocks{j} = chi(chi_index) * N1_blocks{j};
        else
            error('Required chi value (chi_%d) not available', chi_index);
        end
    end
    N2 = vertcat(N2_blocks{:});
    
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
            if L_index > length(chi)
                error('Required chi value (chi_%d) not available for W construction', k+p);
            end
            
            coeff = factorial(k + p - 1)*sqrt(k * p) / (factorial(k) * factorial(p));
            W_block = coeff * chi(k) * R1^(k + p) * L{k + p} ;
            
            % Place the block in the correct position
            row_range = (2*s-1):(2*s);
            col_range = (2*t-1):(2*t);
            W(row_range, col_range) = W_block;
        end
    end
    
    % Calculate Z matrix
    I2 = eye(2);
    Z = J1 - chi(1) * N1 * inv(eye(2*(m)) + W) * N2;
 
    
    % Calculate effective conductivity coefficients
    kappa_eff = coeficientes_efectivos(rho1,rho2,V3,V2,Z);
          
end



%**************************************************************************
% FORMULA DE LOS COEFICIENTES EFECTIVOS (PROBLEMAS Antiplano)
%**************************************************************************
function RR =coeficientes_efectivos(x1,x2,V3,V2,Z)

Delta= (V2*(1-x1)*(x1+x2)+2*V3*x1*(1-x2)) / ...
         ((1+x1)*(x1+x2)*(V2+V3)+(1-x1)*(x1-x2)*(V3));
DtZ=det(Z);
c55=(1-2*(V2+V3)*Delta*(Z(2,2))/(DtZ));
  c45=2*((V2+V3)*Delta*(Z(2,1))/(DtZ));
  c54=2*((V2+V3)*Delta*(Z(1,2))/(DtZ));
c44=(1-2*(V2+V3)*Delta*(Z(1,1))/(DtZ));

      
RR=[c55,c44,c54,c45];
end
