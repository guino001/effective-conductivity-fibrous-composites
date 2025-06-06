clear all
clf

% ii=input('Un numero para el �ngulo: ');
ii=7;
% %%% ejes principales
% %%%%% ejes no principales
load S10_45_75_90_delta1_2;
S1=S10_45_75_90_delta1_2(1:36,ii);
Periodos=cos(theta(ii))+1i*sin(theta(ii)); % w2
w2=Periodos;
w1=1;

delta1=2*ZZ(ii);
delta2=2*ZZ1(ii);   %delta1Leg.*Periodos-2*pi*1i;  

% delta1Leg=2*ZZ(ii);
% delta2Leg=2*ZZ1(ii);   %delta1Leg.*Periodos-2*pi*1i;  
  
%%No =5; % *** Orden del sistema ******************************************
% ************************************************************************
V = abs(w1)*abs(w2)*sin(theta(ii)); % *** Volumen de la celda periodica *******
rm=min(sin(theta(ii)/2),0.5);
percolacion =pi*rm.^2./sin(theta(ii));
angulo=['$\theta=$', num2str(theta(ii)*180/pi,3),'$^o$'];
% 
H1=(conj(delta1)*conj(w2)-conj(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
H2=((delta1)*conj(w2)-(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);

%% Figure 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bi=[150];
xas1=[ 60,100 , 150, 500];
% cf=50; %%% propiedad de conductividad de las fibras
cm=1;
% p=[0:.001:.1, 0.1:.01:1];
% alpha=(1-p)./(p*(.1-1)+1);
%alpha=[0:.01:1];
alpha=[0:.001:.1, 0.1:.01:1];
phi2=[.5,.6, .61:.01:.64,.65]+.25;
%phi2=[.75,.83, .84:.01:.89];
%  phi2=[.1,.2,.25,.3,.4,.5];
phi=phi2;
  K=Bi; % par�metro de imperfeccion
 for mm=1:length(xas1)
  
  cf=xas1(mm);
  cm=1;
%   K=cf/(cf-1); % par�metro de imperfeccion
  for tt=1:length(alpha)
      for   jj=1:length(phi)
r1=sqrt(V*phi(jj)/pi);  % radio de la fibra interna
%       V2masV3=phi(jj);
V2=phi(jj);
xas=cf/cm;
X1=(xas - K*(xas - 1))/(xas + K*(xas + 1));

 JJ=eye(2)+X1*r1^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)];
 
      Vc=alpha(tt)*phi(jj);
      Vf=phi(jj)*(1-alpha(tt))/(1-Vc);
       
%       lamCH=ahm3f_elastico_p13_p23_Paralelogramo_Agosto(cm,cf,cmf,V3,V2,V,S1,H1,H2,4);
%        lamCHLarga(tt,jj,mm)=lamCH(1);
       
      lamCH = ahm_imperfect_Paralelogramo_Cortas_2023_No6(cm,cf,K,V2,V,S1,H1,H2); % Compuesto bifasico 
%       ahm_imperfect_Paralelogramo_Cortas_2023_No4(cm,cf,K,V2,S1,X1,r1,JJ);
%       
            lamInt =ahm_imperfect_Paralelogramo_Cortas_2023_No6(cm,cf,K,Vf,V,S1,H1,H2);
            %RR =AHMelastico_IMPERFECT_sqr_orden6a(cm,cf,phi,K)
            lamRH = ahm_imperfect_Paralelogramo_Cortas_2023_No6(lamInt(1),cf,K,Vc,V,S1,H1,H2);
       lamGain2(tt,jj,mm)=lamRH(1)/lamCH(1); 

  
      end
    end
 end
 
%%%% Figura 5 del paper Julian Ernesto
marca1=[':b';':y';':g';':r';':c';':m';':b';':y';':m'];
marca=['-b';'-y';'-g';'-r';'-c';'-m';'-b';'-y';'-m'];

colores=[[0,0,1];[1,1,0];[0,1,0];[1,0,0];[0,1,1];[1,0,1];[0,0,0];[0.5,0.5,0.5];[0.5,0.5,1];[0.8,0.8,0];[0.1,0.8,0];[.8,0,0];[0,.8,.8]];
xas=xas1;

 
figure(1)
 for mm=1:length(xas)
subplot(2,2,mm)
hold on
for jj=1:length(phi)
     % Seleccionar estilo basado en el �ndice kk
    switch mod(jj,8)
        case 1
            linestyle = '--';  linewidth = 1.5;
        case 2
            linestyle = '-'; linewidth = 1.5;
        case 3
            linestyle = ':';  linewidth = 2.0;
        case 4
            linestyle = '-.'; linewidth = 1.5;
        case 5
            linestyle = '-';  linewidth = 2.5;
        case 6
            linestyle = '--'; linewidth = 2.0;
        case 7
            linestyle = ':';  linewidth = 2.5;
        case 8
            linestyle = '-.';  linewidth = 2.5;
    end
h(jj) =plot(alpha,lamGain2(:,jj,mm),...
                      'LineStyle', linestyle, ...
                       'LineWidth', linewidth, ...
                       'Color', colores(jj,:));

legendEntries{jj} = ['\textbf{$\phi=$}',num2str(phi(jj),3)];
end
plot(alpha,ones(size(alpha)),':k','LineWidth',1,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',4)
% for jj=1:length(phi)                   
% text(.4,lamGain2(136,jj,mm),['\textbf{$\phi=$}',num2str(phi(jj),3)],'fontsize',12,'interpreter','latex')
% end
hold off
% ylim([0.98, 1.08])
xlim([0, 1])
box on
grid on
grid minor
if mm==1 | mm==3
ylabel('$\kappa_{gain}$','fontsize',20,'FontWeight','bold','interpreter','latex')
end
if mm==3 | mm==4
xlabel('aggregation parameter ($\alpha$)','fontsize',16,'FontWeight','bold','interpreter','latex')
end
title( ['$\rho=$',num2str(xas(mm),3),',    ','$Bi=$',num2str(Bi,3),', ',angulo],...
                    'fontsize',12,'interpreter','latex')
                
if mm==1
    % Leyenda
    legend(h, legendEntries, 'FontSize', 12, 'Location', 'best', 'interpreter', 'latex', 'FontWeight', 'bold','NumColumns',7);
end
 end



%%%FUNCIONES
%%% FUNCIONES AUXILIARES

%%% orden del sistema No=1
function RR1 =ahm_imperfect_Paralelogramo_Cortas_2023_No1(cm,cf,K,V2,V,S1,H1,H2)
r1=sqrt(V*V2/pi);  % radio de la fibra interna
xas=cf/cm;
X1=(xas - K*(xas - 1))/(xas + K*(xas + 1));
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
J=eye(2)+X1*r1^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)]; 

S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);S16=S1(16);S18=S1(18);
S20=S1(20);S22=S1(22);S24=S1(24);


X3=(3*xas - K*(xas - 1))/(3*xas + K*(xas + 1));
X5=(5*xas - K*(xas - 1))/(5*xas + K*(xas + 1));
X7=(7*xas - K*(xas - 1))/(7*xas + K*(xas + 1));
X9=(9*xas - K*(xas - 1))/(9*xas + K*(xas + 1));

 
        
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

RR1 =coeficientes_efectivos(cm,V2,X1,Z);
end



%%% orden del sistema No=2
function RR1 =ahm_imperfect_Paralelogramo_Cortas_2023_No2(cm,cf,K,V2,V,S1,H1,H2)
r1=sqrt(V*V2/pi);  % radio de la fibra interna
xas=cf/cm;
X1=(xas - K*(xas - 1))/(xas + K*(xas + 1));
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
J=eye(2)+X1*r1^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)]; 


S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);S16=S1(16);S18=S1(18);
S20=S1(20);S22=S1(22);S24=S1(24);

X1=(xas - K*(xas - 1))/(xas + K*(xas + 1));
X3=(3*xas - K*(xas - 1))/(3*xas + K*(xas + 1));
X5=(5*xas - K*(xas - 1))/(5*xas + K*(xas + 1));
X7=(7*xas - K*(xas - 1))/(7*xas + K*(xas + 1));
X9=(9*xas - K*(xas - 1))/(9*xas + K*(xas + 1));
 
X11=(11*xas - K*(xas - 1))/(11*xas + K*(xas + 1));
X13=(13*xas - K*(xas - 1))/(13*xas + K*(xas + 1));
        
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

RR2 =coeficientes_efectivos(cm,V2,X1,Z);

end


%%% orden del sistema No=3
function RR1 =ahm_imperfect_Paralelogramo_Cortas_2023_No3(cm,cf,K,V2,V,S1,H1,H2)
r1=sqrt(V*V2/pi);  % radio de la fibra interna
xas=cf/cm;
X1=(xas - K*(xas - 1))/(xas + K*(xas + 1));
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
J=eye(2)+X1*r1^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)]; 


S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);S16=S1(16);S18=S1(18);
S20=S1(20);S22=S1(22);S24=S1(24);

X1=(xas - K*(xas - 1))/(xas + K*(xas + 1));
X3=(3*xas - K*(xas - 1))/(3*xas + K*(xas + 1));
X5=(5*xas - K*(xas - 1))/(5*xas + K*(xas + 1));
X7=(7*xas - K*(xas - 1))/(7*xas + K*(xas + 1));
X9=(9*xas - K*(xas - 1))/(9*xas + K*(xas + 1));
 
X11=(11*xas - K*(xas - 1))/(11*xas + K*(xas + 1));
X13=(13*xas - K*(xas - 1))/(13*xas + K*(xas + 1));
        
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

RR3 =coeficientes_efectivos(cm,V2,X1,Z);

end


%%% orden del sistema No=4
function RR4 =ahm_imperfect_Paralelogramo_Cortas_2023_No4(cm,cf,K,V2,V,S1,H1,H2)
r1=sqrt(V*V2/pi);  % radio de la fibra interna
xas=cf/cm;
X1=(xas - K*(xas - 1))/(xas + K*(xas + 1));
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
J=eye(2)+X1*r1^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)]; 

S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);S16=S1(16);S18=S1(18);
S20=S1(20);S22=S1(22);S24=S1(24);

X1=(xas - K*(xas - 1))/(xas + K*(xas + 1));
X3=(3*xas - K*(xas - 1))/(3*xas + K*(xas + 1));
X5=(5*xas - K*(xas - 1))/(5*xas + K*(xas + 1));
X7=(7*xas - K*(xas - 1))/(7*xas + K*(xas + 1));
X9=(9*xas - K*(xas - 1))/(9*xas + K*(xas + 1));
 
X11=(11*xas - K*(xas - 1))/(11*xas + K*(xas + 1));
X13=(13*xas - K*(xas - 1))/(13*xas + K*(xas + 1));
        
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

N11=sqrt(3)*r1^4*L4; N12=sqrt(5)*r1^6*L6; N13=sqrt(7)*r1^8*L8; N14=3*r1^10*L10;

N1=[N11, N12, N13,N14];

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

RR4 =coeficientes_efectivos(cm,V2,X1,Z);

end


%%% orden del sistema No=5
function RR5 =ahm_imperfect_Paralelogramo_Cortas_2023_No5(cm,cf,K,V2,V,S1,H1,H2)
r1=sqrt(V*V2/pi);  % radio de la fibra interna
xas=cf/cm;
X1=(xas - K*(xas - 1))/(xas + K*(xas + 1));
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
J=eye(2)+X1*r1^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)]; 

S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);S16=S1(16);S18=S1(18);
S20=S1(20);S22=S1(22);S24=S1(24);

X1=(xas - K*(xas - 1))/(xas + K*(xas + 1));
X3=(3*xas - K*(xas - 1))/(3*xas + K*(xas + 1));
X5=(5*xas - K*(xas - 1))/(5*xas + K*(xas + 1));
X7=(7*xas - K*(xas - 1))/(7*xas + K*(xas + 1));
X9=(9*xas - K*(xas - 1))/(9*xas + K*(xas + 1));
 
X11=(11*xas - K*(xas - 1))/(11*xas + K*(xas + 1));
X13=(13*xas - K*(xas - 1))/(13*xas + K*(xas + 1));
        
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

RR5 =coeficientes_efectivos(cm,V2,X1,Z);

end


%%% orden del sistema No=6
function RR6 =ahm_imperfect_Paralelogramo_Cortas_2023_No6(cm,cf,K,V2,V,S1,H1,H2)
r1=sqrt(V*V2/pi);  % radio de la fibra interna
xas=cf/cm;
X1=(xas - K*(xas - 1))/(xas + K*(xas + 1));
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
J=eye(2)+X1*r1^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)]; 

S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);S16=S1(16);S18=S1(18);
S20=S1(20);S22=S1(22);S24=S1(24);S26=S1(26);

X1=(xas - K*(xas - 1))/(xas + K*(xas + 1));
X3=(3*xas - K*(xas - 1))/(3*xas + K*(xas + 1));
X5=(5*xas - K*(xas - 1))/(5*xas + K*(xas + 1));
X7=(7*xas - K*(xas - 1))/(7*xas + K*(xas + 1));
X9=(9*xas - K*(xas - 1))/(9*xas + K*(xas + 1));
 
X11=(11*xas - K*(xas - 1))/(11*xas + K*(xas + 1));
X13=(13*xas - K*(xas - 1))/(13*xas + K*(xas + 1));
        
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
W41=55*sqrt(3)*X9*r1^12*L12;  W42=143*sqrt(45)*X9*r1^14*L14;   W43=2145*sqrt(7)*X9*r1^16*L16;   W44=24310*X9*r1^18*L18;           W45=25194*sqrt(11)*X9*r1^20*L20; 
W51=26*sqrt(33)*X11*r1^14*L14; W52=273*sqrt(55)*X11*r1^16*L16; W53=1768*sqrt(77)*X11*r1^18*L18; W54=25194*sqrt(11)*X11*r1^20*L20; W55=352716*X11*r1^22*L22;
W61=35*sqrt(39)*X13*r1^16*L16; W62=467*sqrt(65)*X13*r1^18*L18; W63=3876*sqrt(91)*X13*r1^20*L20; W64=22610*sqrt(117)*X13*r1^22*L22; W65=104006*sqrt(143)*X13*r1^24*L24;
W16=35*sqrt(39)*X3*r1^16*L16; W26=467*sqrt(65)*X5*r1^18*L18;  W36=3876*sqrt(91)*X7*r1^20*L20; W46=22610*sqrt(117)*X9*r1^22*L22; W56=104006*sqrt(143)*X11*r1^24*L24;  
W66=5200300*X13*r1^26*L26;


W=[W11, W12, W13, W14, W15, W16; W21, W22, W23, W24, W25, W26; W31, W32, W33, W34, W35, W36;...
   W41, W42, W43, W44, W45, W46; W51, W52, W53, W54, W55, W56; W61, W62, W63, W64, W65, W66];

AP2=eye(12)+W;
PP=-X1*N1*inv(AP2)*N2;  
JK=J+PP; % Matriz Z
% disp('corto');
Z=JK;

RR6 =coeficientes_efectivos(cm,V2,X1,Z);

end


function [kappa_eff] = conductividad_efectiva_imperfect_orden_m(cm,cf,K,V2,V,Sk,H1,H2,m)

R1=sqrt(V*V2/pi);  % radio de la fibra interna
rho=cf/cm;
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
 

    % Calculate chi values up to required order
    max_p = 6*(m) + 1; % We need up to chi_{2(m)+1}
    chi = zeros(1, max_p);
    for p = 1:2:max_p
        numerator = p*rho-K*(rho-1);
        denominator = p*rho+K*(rho+1);
        chi(p) = numerator / denominator;
    end
  J1=eye(2)+chi(1)*R1^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)];     
  
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
    kappa_eff = coeficientes_efectivos(cm,V2,chi(1),Z);
end

%**************************************************************************
% FORMULA DE LOS COEFICIENTES EFECTIVOS (PROBLEMAS Antiplano)
%**************************************************************************
function RR =coeficientes_efectivos(cm,V2,X1,Z)
DtZ=det(Z);
c55=cm*(1-2*V2*X1*Z(2,2)/DtZ);
  c45=cm*(2*V2*X1*Z(2,1)/DtZ);
  c54=cm*(2*V2*X1*Z(1,2)/DtZ);
c44=cm*(1-2*V2*X1*Z(1,1)/DtZ);
      
RR=[c55,c44,c54,c45];
end

 
 
 function RR=ahm_imperfect_Paralelogramo(cm,cf,K,V2,V,S1,H1,H2,No)  
 
r = sqrt(V*V2/pi);  % *** Radio de la fibra *******************************
xas=cf/cm;
S=S1;
%Matriz Wnk
for n=1:2:30
   for k=1:2:30

Wnk(n,k)=(S(n+k))*(r^(n+k))*((factorial(n+k-1)/(factorial(n-1)*factorial(k-1)))/...
               sqrt(n*k));

   end
end
RWnk=real(Wnk);
IWnk=imag(Wnk);
di=eye(70); 

% Para hallar los Bp

for s=1:2:30,

       Bp(s)=((1-xas)*K+s*xas)/((xas+1)*K+ s*xas);
       
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