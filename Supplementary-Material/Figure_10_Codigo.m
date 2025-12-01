clearvars
clear all

 %%%%%%%%%%%%%%%%%%%%%
% 
load S10_45_75_90_delta1_2; 

cm=1;
ii=11;
% ************************************************************************
S1=S10_45_75_90_delta1_2(:,ii);
 delta1=2*ZZ(ii);
 Periodos2=(cos(theta(ii))+1i*sin(theta(ii)));
 w1=1;
 w2=Periodos2;
rm=min(sin(theta(ii)/2),0.5);
percolacion =pi*rm.^2./sin(theta(ii));

  
delta1=delta1;
delta2=delta1.*w2-2*pi*1i; 

V = abs(w1)*abs(w2)*sin(theta(ii)); % *** Volumen de la celda periodica *******
% 
%H1=(conj(delta1)*conj(w2)-conj(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
H2=((delta1)*conj(w2)-(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
H1=(pi)/imag(w1*w2);
% h11=real(H1);
% h21=imag(H1);
h12=real(H2);
h22=imag(H2);
phi=[.785]; 
% *************************************************************************
alpha=[0.01,.1,.2,.3:.1:.9,.95];
xas=[1:1:100];
cf=xas;
X=(1-xas)./(1+xas);
for tt=1:length(alpha)
for jj=1:length(xas)
 Vc=alpha(tt)*phi;
 Vf=phi*(1-alpha(tt))/(1-Vc);
 %% AHM n=0;
% lamInt=ahm2f_Paralelogramo_Cortas_2023_No0(cm,cf(jj),Vf,V,S1,H1,H2); % VanFo Fi
% CENo(jj,:,tt)=ahm2f_Paralelogramo_Cortas_2023_No0(lamInt(1),cf(jj),Vc,V,S1,H1,H2); % VanFo Fi
%% AHM n=1;
lamInt =ahm2f_Paralelogramo_Cortas_2023_No1(cm,cf(jj),Vf,V,S1,H1,H2);
 CENo(jj,:,tt) =ahm2f_Paralelogramo_Cortas_2023_No1(lamInt(1),cf(jj),Vc,V,S1,H1,H2);


%% AHM n=2;
% lamInt =ahm2f_Paralelogramo_Cortas_2023_No2(cm,cf,Vf,V,S1,H1,H2);
% CENo(jj,:,tt) =ahm2f_Paralelogramo_Cortas_2023_No2(lamInt(1),cf(jj),Vc,V,S1,H1,H2);

%  CENo3(j,:,tt) =ahm2f_Paralelogramo_Cortas_2023_No3(cm,cf,V2,V,S1); 
%  CENo4(j,:,tt) =ahm2f_Paralelogramo_Cortas_2023_No4(cm,cf,V2,V,S1); 

%  lamInt =ahm2f_Paralelogramo_Cortas_2023_No5(cm,cf(jj),Vf,V,S1,H1,H2);
%  CENo(jj,:,tt) =ahm2f_Paralelogramo_Cortas_2023_No5(lamInt(1),cf(jj),Vc,V,S1,H1,H2);

%  lamInt =conductividad_efectiva_orden_m(cm,cf(jj),Vf,V,S1,H1,H2,30);
%  CENo(jj,:,tt) =conductividad_efectiva_orden_m(lamInt(1),cf(jj),Vc,V,S1,H1,H2,10);
end                                 
end
    
angulo(tt,:)=['$\theta=$', num2str(theta(ii)*180/pi,3),'$^o$'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


matriz=[ '$\kappa_{1}=$',num2str(cm,5)];
fibra=[ '$\kappa_{2}=$',num2str(cf,5)];
Modulo=[ '|w2|= ',num2str(abs(w2),5)];
roo=[ '$\rho=$',num2str(cf/cm,5)];

%textoY=['$\hat{\kappa}_{11}/\kappa_1$';'$\hat{\kappa}_{22}/\kappa_1$';'$\hat{\kappa}_{12}/\kappa_1$'];
disp(Modulo)
disp(angulo(tt,:))
disp(matriz)
disp(fibra)
%Tabla=zeros(length(phi),4,length(alpha));
Tabla(:,:)=[CENo(:,:,1)];


 %%% Cotas de Hashin
 f=phi;
  HS1=cm.*(cf.*(1+f)+cm.*(1-f))./(cf.*(1-f)+cm.*(1+f));
  HS2=cf.*(cf.*f +cm.*(2-f))./(cf.*(2-f)+cm.*f);
  
  if cf(1)<=1 
    HU=HS2;
    HL=HS1;
else
    HU=HS1;
    HL=HS2;  
end
 %*******************************************************************************
% FORMULA DE LAS COTAS DE BRUNO, 1991
%*********************************************************************
z=cf;p1=phi;d=2;
p2=1-p1; % (matrix volume fraction)
q=2*sqrt(p1.*V/pi);
sm=0.5*(1-q.^2);
SM=0.5*(1+q.^2);
delta=SM-sm;
V=1-phi./((1./(1-z))-((1-phi)./2)); %formula (80) p. 369
W=(1-(1-z).*SM)./(1-(1-z).*sm).*(1+(1-p1./delta).^2./((1-p1./delta).*...
 (1-(1-z).*SM)./((1-z).*delta)+(p1./delta.^2).*(p2./d-sm)));       
Gaf=29.959016;
Gam=1.2777778;

if z(1)<=1 
    BU=W;
    BL=V;
else
    BU=V;
    BL=W;   
end

% Datos experimentales pag. 60 de Aboudi

XGA=[0.48 0.52 0.546 0.56 0.6 0.656 0.672 0.62 0.65 0.66];

YGA=[2.4615 3.0769 3.1538 3.3077 3.2308 3.6308 4.0461 ...
     4.3846 4.6923 4.8461];
% UGa=Gaf+(1-p1(jj))/(1/(Gam-Gaf)+p1(jj)/(2*Gaf));
	       
%*********************************************************
 f=phi; rho=-X;
%MIT = MITYUSHEV_sqr_100(f, rho);
MIT =point7_upper_bound(cf, cm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

colores = [[0,0,1];[.4,.1,.4];[0,1,0];[1,0,0];[0,1,1];...
    [1,0,.8];[.8,0.5,.8];[0.5,0.5,0.5];[0.5,0.5,1];[0.8,0,0.7];[0.1,0.8,0];[.8,0,0];[0,.7,.5];[.6,.3,0];[0.2,1,.7];[4,0,0]];
marca = ['ko-'; 'ko-';'ko-';'ko-'; 'ko-';'ko-';'ko-';'ks-'];
letras = ['a)';'b)';'c)';'d)';'e)';'f)'];

figure(1)
subplot(1,2,1)
hold on

% Crear handles para la leyenda
h = [];

% Primer plot - HU
h1 = plot(xas, HU, 'k-', 'LineWidth', 2, 'MarkerEdgeColor', 'k', ...
          'MarkerFaceColor', 'w', 'MarkerSize', 6);
h = [h, h1];
% Plot MIT
h2 = plot(xas, MIT, 'b-.', 'LineWidth', 2, 'MarkerEdgeColor', 'r', ...
          'MarkerFaceColor', 'w', 'MarkerSize', 6);
h = [h, h2];
% Plots para diferentes valores de alpha

h3 = plot(xas, BU, 'r-.', 'LineWidth', 2.5, 'MarkerEdgeColor', 'r', ...
          'MarkerFaceColor', 'w', 'MarkerSize', 6);
h = [h, h3];
for tt = 1:length(alpha)
    h_temp = plot(xas, CENo(:,1,tt),'LineWidth', 2, 'Color', colores(tt,:));
    h = [h, h_temp];
end
h3 = plot(xas, BU, 'r-.', 'LineWidth', 2.5, 'MarkerEdgeColor', 'r', ...
          'MarkerFaceColor', 'w', 'MarkerSize', 6);


h4 = plot(xas, BL, 'r--', 'LineWidth', 2.5, 'MarkerEdgeColor', 'r', ...
          'MarkerFaceColor', 'w', 'MarkerSize', 6);
h = [h, h4];

% Plot HL
h5 = plot(xas, HL, 'k-.', 'LineWidth', 2, 'MarkerEdgeColor', 'r', ...
          'MarkerFaceColor', 'w', 'MarkerSize', 6);
h = [h, h5];
hold off
text(xas(3),HU(3),sprintf('$V_2$ = %.3f',phi),'fontsize', 14, 'interpreter', 'latex')
% Crear las etiquetas para la leyenda
legend_labels = {'Hashin Shtrikman upper bound'};
legend_labels{end+1} = '$\sigma_{7}^U$';
legend_labels{end+1} = 'O.P. Bruno upper bound';
% Agregar etiquetas para los plots de alpha
for tt = 1:length(alpha)
    legend_labels{end+1} = sprintf('$\\alpha = %.2f$', alpha(tt));
end

% Agregar las dem·s etiquetas

legend_labels{end+1} = 'O.P. Bruno lower bound';
legend_labels{end+1} = 'Hashin Shtrikman lower bound';
% Configurar la leyenda
legend(h, legend_labels, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 10);

xlabel('$\kappa_2/\kappa_1$', 'fontsize', 18, 'interpreter', 'latex')
xlim([min(xas), max(xas)])
a=round((max(xas)-min(xas))/11);
%xticks([min(xas):a:max(xas)])
xticks([1,10:10:100])
ylabel('$\hat{\kappa}/\kappa_1$', 'fontsize', 18, 'interpreter', 'latex')
grid on
grid minor
box on


subplot(1,2,2)
hold on

% Crear handles para la leyenda
h = [];

% Primer plot - HU
h1 = plot(xas, HU, 'k-', 'LineWidth', 2, 'MarkerEdgeColor', 'k', ...
          'MarkerFaceColor', 'w', 'MarkerSize', 6);
h = [h, h1];
% Plot MIT
h2 = plot(xas, MIT, 'b-.', 'LineWidth', 2, 'MarkerEdgeColor', 'r', ...
          'MarkerFaceColor', 'w', 'MarkerSize', 6);
h = [h, h2];
% Plots para diferentes valores de alpha

h3 = plot(xas, BU, 'r-.', 'LineWidth', 2.5, 'MarkerEdgeColor', 'r', ...
          'MarkerFaceColor', 'w', 'MarkerSize', 6);
h = [h, h3];
for tt = 1:length(alpha)
    h_temp = plot(xas, CENo(:,1,tt),'LineWidth', 2, 'Color', colores(tt,:));
    h = [h, h_temp];
end
h3 = plot(xas, BU, 'r-.', 'LineWidth', 2.5, 'MarkerEdgeColor', 'r', ...
          'MarkerFaceColor', 'w', 'MarkerSize', 6);


h4 = plot(xas, BL, 'r--', 'LineWidth', 2.5, 'MarkerEdgeColor', 'r', ...
          'MarkerFaceColor', 'w', 'MarkerSize', 6);
h = [h, h4];

% Plot HL
h5 = plot(xas, HL, 'k-.', 'LineWidth', 2, 'MarkerEdgeColor', 'r', ...
          'MarkerFaceColor', 'w', 'MarkerSize', 6);
h = [h, h5];
hold off
text(xas(3),HU(3),sprintf('$V_2$ = %.3f',phi),'fontsize', 14, 'interpreter', 'latex')
% Crear las etiquetas para la leyenda
legend_labels = {'Hashin Shtrikman upper bound'};
legend_labels{end+1} = '$\sigma_{7}^U$';
legend_labels{end+1} = 'O.P. Bruno upper bound';
% Agregar etiquetas para los plots de alpha
for tt = 1:length(alpha)
    legend_labels{end+1} = sprintf('$\\alpha = %.2f$', alpha(tt));
end

% Agregar las dem·s etiquetas

legend_labels{end+1} = 'O.P. Bruno lower bound';
legend_labels{end+1} = 'Hashin Shtrikman lower bound';
% Configurar la leyenda
legend(h, legend_labels, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 10);

xlabel('$\kappa_2/\kappa_1$', 'fontsize', 18, 'interpreter', 'latex')
xlim([1, 10])
xticks([1:1:10])
ylim([1, 7])
ylabel('$\hat{\kappa}/\kappa_1$', 'fontsize', 18, 'interpreter', 'latex')
grid on
grid minor
box on

% Funciones auxiliares
% Funciones auxiliares
%% para VanFo Fi
function RR1=ahm2f_Paralelogramo_Cortas_2023_No0(cm,cf,V2,V,S1,H1,H2)
 r1=sqrt(V*V2/pi);  % radio de la fibra interna
x=cf/cm;
X=(1-x)/(1+x);
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
J=eye(2)+X*r1^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)];
RR1 =coeficientes_efectivos2f(cm,V2,X,J); % VanFo Fi
%%% orden del sistema No=1
end
%%% orden del sistema No=1
function RR1=ahm2f_Paralelogramo_Cortas_2023_No1(cm,cf,V2,V,S1,H1,H2)
r1=sqrt(V*V2/pi);  % radio de la fibra interna
x=cf/cm;
X=(1-x)/(1+x);
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
J=eye(2)+X*r1^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)];

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
function RR2 =ahm2f_Paralelogramo_Cortas_2023_No2(cm,cf,V2,V,S1,H1,H2)
r1=sqrt(V*V2/pi);  % radio de la fibra interna
x=cf/cm;
X=(1-x)/(1+x);
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
J=eye(2)+X*r1^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)];
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
function RR3 =ahm2f_Paralelogramo_Cortas_2023_No3(cm,cf,V2,V,S1,H1,H2)
r1=sqrt(V*V2/pi);  % radio de la fibra interna
x=cf/cm;
X=(1-x)/(1+x);
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
J=eye(2)+X*r1^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)];
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


function RR3 =ahm2f_Paralelogramo_Cortas_2023_No4(cm,cf,V2,V,S1,H1,H2)
r1=sqrt(V*V2/pi);  % radio de la fibra interna
x=cf/cm;
X=(1-x)/(1+x);
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
J=eye(2)+X*r1^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)];

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
W31=12*sqrt(21)*r1^10*L10;  W32=66*sqrt(35)*r1^12*L12;   W33=1716*r1^14*L14;           W34=2145*sqrt(7)*r1^16*L16;     
W41=55*sqrt(3)*r1^12*L12;  W42=429*sqrt(5)*r1^14*L14;   W43=2145*sqrt(7)*r1^16*L16;   W44=24310*r1^18*L18;  

W=[W11, W12, W13, W14; W21, W22, W23, W24;W31, W32, W33, W34; W41, W42, W43, W44];
    
AP2=eye(8)+X*W;
PP=-X^2*N1*inv(AP2)*N2;  
JK=J+PP; % Matriz Z
Z=JK;

RR3 =coeficientes_efectivos2f(cm,V2,X,Z);

end

%%% orden del sistema No=5
function RR3 =ahm2f_Paralelogramo_Cortas_2023_No5(cm,cf,V2,V,S1,H1,H2)
r1=sqrt(V*V2/pi);  % radio de la fibra interna
x=cf/cm;
X=(1-x)/(1+x);
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
J=eye(2)+X*r1^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)];

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
function RR6 =ahm2f_Paralelogramo_Cortas_2023_No6(cm,cf,V2,V,S1)
r1=sqrt(V*V2/pi);  % radio de la fibra interna
x=cf/cm;
X=(1-x)/(1+x);
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
J=eye(2)+X*r1^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)];

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


function [kappa_eff] = conductividad_efectiva_orden_m(cm,cf,V2,V,Sk,H1,H2,m)
R=sqrt(V*V2/pi);  % radio de la fibra interna
x=cf/cm;
X=(1-x)/(1+x);
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
J=eye(2)+X*R^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)];
                    
                    
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

function sigma_ratio = MITYUSHEV_sqr_100(f, rho)
% CALCULATE_SIGMA_RATIO Calcula sigma_100(f,rho)/sigma_2 seg√∫n la f√≥rmula proporcionada
%
% Inputs:
%   f   - par√°metro libre f (puede ser escalar o vector)
%   rho - par√°metro libre rho (puede ser escalar o vector)
%
% Output:
%   sigma_ratio - valor de sigma_100(f,rho)/sigma_2

    % Inicializar el resultado
    sigma_ratio = zeros(size(f));
    
    % TÈrminos de la serie
    for i = 1:numel(rho)
        current_f = f;
        current_rho = rho(i);
        
        % TÈrminos principales
        term = 1 + 2*current_f*current_rho + 2*current_f^2*current_rho^2;
        
        % TÈrminos con rho^3
        term = term + 2.30974*current_f^3*current_rho^3;
        term = term + 0.358151*current_f^4*current_rho^3;
        term = term + 0.440133*current_f^5*current_rho^3;
        term = term + 0.283968*current_f^6*current_rho^3;
        term = term + 0.17321*current_f^7*current_rho^3;
        term = term + 0.0763706*current_f^8*current_rho^3;
        term = term + 0.0349588*current_f^9*current_rho^3;
        term = term + 0.0192229*current_f^10*current_rho^3;
        term = term + 0.0100547*current_f^11*current_rho^3;
        term = term + 0.00455705*current_f^12*current_rho^3;
        
        % TÈrminos con rho^4
        term = term + 2.63336*current_f^4*current_rho^4;
        term = term + 0.764796*current_f^5*current_rho^4;
        term = term + 0.828513*current_f^6*current_rho^4;
        term = term + 0.591477*current_f^7*current_rho^4;
        term = term + 0.357103*current_f^8*current_rho^4;
        term = term + 0.281438*current_f^9*current_rho^4;
        term = term + 0.255881*current_f^10*current_rho^4;
        term = term + 0.175048*current_f^11*current_rho^4;
        term = term + 0.0607134*current_f^12*current_rho^4;
        
        % TÈrminos con rho^5
        term = term + 3.20538*current_f^5*current_rho^5;
        term = term + 1.42227*current_f^6*current_rho^5;
        term = term + 1.71443*current_f^7*current_rho^5;
        term = term + 1.69485*current_f^8*current_rho^5;
        term = term + 1.62367*current_f^9*current_rho^5;
        term = term + 2.09162*current_f^10*current_rho^5;
        term = term + 2.7258*current_f^11*current_rho^5;
        term = term + 3.72375*current_f^12*current_rho^5;
        
        % TÈrminos con rho^6
        term = term + 3.87034*current_f^6*current_rho^6;
        term = term + 2.25248*current_f^7*current_rho^6;
        term = term + 2.86381*current_f^8*current_rho^6;
        term = term + 3.14345*current_f^9*current_rho^6;
        term = term + 3.1104*current_f^10*current_rho^6;
        term = term + 4.40338*current_f^11*current_rho^6;
        term = term + 5.4525*current_f^12*current_rho^6;
        
        % TÈrminos con rho^7
        term = term + 4.87402*current_f^7*current_rho^7;
        term = term + 3.5852*current_f^8*current_rho^7;
        term = term + 4.91752*current_f^9*current_rho^7;
        term = term + 5.87811*current_f^10*current_rho^7;
        term = term + 6.99559*current_f^11*current_rho^7;
        term = term + 10.3469*current_f^12*current_rho^7;
        
        % TÈrminos con rho^8
        term = term + 6.07698*current_f^8*current_rho^8;
        term = term + 5.34872*current_f^9*current_rho^8;
        term = term + 7.69997*current_f^10*current_rho^8;
        term = term + 9.75023*current_f^11*current_rho^8;
        term = term + 12.5134*current_f^12*current_rho^8;
        
        % TÈrminos con rho^9
        term = term + 7.80964*current_f^9*current_rho^9;
        term = term + 8.07344*current_f^10*current_rho^9;
        term = term + 12.3949*current_f^11*current_rho^9;
        term = term + 16.5787*current_f^12*current_rho^9;
        
        % TÈrminos con rho^10
        term = term + 9.93356*current_f^10*current_rho^10;
        term = term + 11.7505*current_f^11*current_rho^10;
        term = term + 18.8899*current_f^12*current_rho^10;
        
        % TÈrminos con rho^11
        term = term + 12.9369*current_f^11*current_rho^11;
        term = term + 17.2854*current_f^12*current_rho^11;
        
        % TÈrminos con rho^12
        term = term + 16.6825*current_f^12*current_rho^12;
        
        sigma_ratio(i) = term;
    end
end

function sigma_100 = point7_upper_bound(sigma1, sigma2)
   
    %%% (cf=sigma1, cm= sigma2);
    
    % TÈrminos DEL NUMERADOR (4 TÈrminos)
    term_n1 = 3.634660* sigma1.^3;
    term_n2 = 13.7384 * sigma1.^2 * sigma2;  %
    term_n3 = 7.83067* sigma1 .* sigma2.^2;
    term_n4 = 0.924842 * sigma2.^3;
    
    numerador = term_n1 + term_n2 + term_n3 + term_n4;
    
    % TÈrminos DEL DENOMINADOR (5 TÈrminos)  
    term_d1 = sigma1.^4;  
    term_d2 = 8.22711 * sigma1.^3 * sigma2; 
    term_d3 = 13.5477 * sigma1.^2* sigma2.^2;    
    term_d4 = 3.35568* sigma1* sigma2.^3;
    term_d5 = -0.0019548* sigma2.^4;
    
    
    denominador = term_d1 + term_d2 + term_d3 + term_d4 + term_d5;
    
    % EVALUAR LA F√ìRMULA COMPLETA
    sigma_100 = sigma1 .* (numerador ./ denominador);
end