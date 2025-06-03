clear all
clf

% ii=input('Un numero para el ángulo: ');
ii=7;
% %%% ejes principales
% load S10_30_45_75_90EP
% Sk=S10_30_45_75_90EP(:,ii);
% Periodos=cos(theta(ii)/2)-1i*sin(theta(ii)/2); % w2
% w2=Periodos;
% w1=conj(w2);


% %%%%% ejes no principales
load S10_45_75_90_delta1_2;
S1=S10_45_75_90_delta1_2(:,ii);
Periodos=cos(theta(ii))+1i*sin(theta(ii)); % w2
w2=Periodos;
w1=1;

delta1=2*ZZ(ii);
delta2=2*ZZ1(ii);   %delta1Leg.*Periodos-2*pi*1i;  
angulo=[', ','$\theta=$',num2str(theta(ii)*180/pi,3),'$^o$'];
% delta1Leg=2*ZZ(ii);
% delta2Leg=2*ZZ1(ii);   %delta1Leg.*Periodos-2*pi*1i;  
  
No =7; % *** Orden del sistema ******************************************
V = abs(w1)*abs(w2)*sin(theta(ii)); % *** Volumen de la celda periodica *******
% *** expresiones de H1 y H2 en el problema I23 ********
H1=(conj(delta1)*conj(w2)-conj(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
H2=((delta1)*conj(w2)-(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
% *************************************************************************
Lh=[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)];
Bi=[15];
xas=[ 60,100 , 150, 500];
% cf=50; %%% propiedad de conductividad de las fibras
cm=1;
% p=[0:.001:.1, 0.1:.01:1];
% alpha=(1-p)./(p*(.1-1)+1);
%alpha=[0:.01:1];
alpha=[0:.001:.1, 0.1:.01:1];
phi2=[.5,.6, .61:.01:.64,.65]+.25; %%V3
V2=0.001;
x1=5;
 for mm=1:length(xas)
  phi=phi2;
  x2=xas(mm);%%% x2=cf/cm pero cm=1
  for ii=1:length(alpha)
      for   jj=1:length(phi)
        r2=sqrt(V*phi(jj)/pi);  % radio de la fibra interna
        V3=phi(jj);
      lamCH = ahm3f_Paralelogramo_Cortas_2025_No5(cm,x1,x2,V3,V2,V,S1,Lh); % Compuesto bifasico
      Vc=alpha(ii)*phi(jj);
      Vf=phi(jj)*(1-alpha(ii))/(1-Vc);
      lamInt =ahm3f_Paralelogramo_Cortas_2025_No5(cm,x1,x2,Vf,V2,V,S1,Lh);
      lamRH = ahm3f_Paralelogramo_Cortas_2025_No5(lamInt(1),x1,x2,Vc,V2,V,S1,Lh);
      lamInt1(ii,jj,mm)=lamInt(1);
      lamRH1(ii,jj,mm)=lamRH(1);
      lamGain2(ii,jj,mm)=lamRH(1)/lamCH(1);
  
  
      end
    end
 end
 
%%%% Figura 5 del paper Julian Ernesto
marca1=[':b';':y';':g';':r';':c';':m';':b';':y';':m'];
marca=['-b';'-y';'-g';'-r';'-c';'-m';'-b';'-y';'-m'];

colores=[[0,0,1];[1,1,0];[0,1,0];[1,0,0];[0,1,1];[1,0,1];[0,0,0];[0.5,0.5,0.5];[0.5,0.5,1];[0.8,0.8,0];[0.1,0.8,0];[.8,0,0];[0,.8,.8]];


figure(1)
 for mm=1:length(xas)
subplot(2,2,mm)
hold on
for jj=1:length(phi)
     % Seleccionar estilo basado en el índice kk
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

legendEntries{jj} = ['\textbf{$V_3=$}',num2str(phi(jj),3)];
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
title( ['$\rho_1=$',num2str(x1,3),',    ','$\rho_2=$',num2str(xas(mm),3), ',    ','$V_2=$',num2str(V2,3),angulo],...
        'fontsize',12,'interpreter','latex')
                
if mm==1
    % Leyenda
    legend(h, legendEntries, 'FontSize', 12, 'Location', 'best', 'interpreter', 'latex', 'FontWeight', 'bold','NumColumns',7);
end
 end
%% END figure 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%% Funciones auxiliares 
%%% orden del sistema No=1
function RR1 =ahm3f_Paralelogramo_Cortas_2025_No1(x1,x2,V3,V2,S1,X1,r1,J)

S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);S16=S1(16);S18=S1(18);
S20=S1(20);S22=S1(22);S24=S1(24);

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
function RR2 =ahm3f_Paralelogramo_Cortas_2025_No2(x1,x2,V3,V2,S1,X1,r1,J)

S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);S16=S1(16);S18=S1(18);
S20=S1(20);S22=S1(22);S24=S1(24);

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
function RR3 =ahm3f_Paralelogramo_Cortas_2025_No3(x1,x2,V3,V2,S1,X1,r1,J)

S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);S16=S1(16);S18=S1(18);
S20=S1(20);S22=S1(22);S24=S1(24);

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
function RR4 =ahm3f_Paralelogramo_Cortas_2025_No4(x1,x2,V3,V2,S1,X1,r1,J)

S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);S16=S1(16);S18=S1(18);
S20=S1(20);S22=S1(22);S24=S1(24);

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

N11=sqrt(3)*r1^4*L4; N12=sqrt(5)*r1^6*L6; N13=sqrt(7)*r1^8*L8; N14=3*r1^10*L10;

N1=[N11, N12, N13,N14];

N2=[X3*N11; X5*N12; X7*N13; X9*N14];

W11=10*X3*r1^6*L6;             W12=7*sqrt(15)*X3*r1^8*L8;      W13=12*sqrt(21)*X3*r1^10*L10;    W14=55*sqrt(27)*X3*r1^12*L12;     
W21=7*sqrt(15)*X5*r1^8*L8;     W22=126*X5*r1^10*L10;           W23=66*sqrt(35)*X5*r1^12*L12;    W24=429*sqrt(5)*X5*r1^14*L14;     
W31=12*sqrt(21)*X7*r1^10*L10;  W32=66*sqrt(35)*X7*r1^12*L12;   W33=1716*X7*r1^14*L14;           W34=2145*sqrt(7)*X7*r1^16*L16;     
W41=55*sqrt(27)*X9*r1^12*L12;  W42=429*sqrt(5)*X9*r1^14*L14;   W43=2145*sqrt(7)*X9*r1^16*L16;   W44=24310*X9*r1^18*L18;             

W=[W11, W12, W13, W14; W21, W22, W23, W24;W31, W32, W33, W34; W41, W42, W43, W44];

AP2=eye(8)+W;
PP=-X1*N1*inv(AP2)*N2;  
JK=J+PP; % Matriz Z
% disp('corto');
Z=JK;

RR4 =coeficientes_efectivos(x1,x2,V3,V2,Z);

end


%%% orden del sistema No=5
function RR5 =ahm3f_Paralelogramo_Cortas_2025_No5(cm,Kt,x2,V3,V2,V,S1,Lh)
r1 = sqrt(V*(V2+V3)/pi);  % *** Radio de la fibra más la mesophase 
x1=Kt/r1; %cmf= x1=p_I/p1
X1=((1-x1)*(x1+x2)*(V2+V3)+(1+x1)*(x1-x2)*(V3))/...
             ((1+x1)*(x1+x2)*(V2+V3)+(1-x1)*(x1-x2)*(V3));   
J=eye(2)+X1*r1^2*Lh;
S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);S16=S1(16);S18=S1(18);
S20=S1(20);S22=S1(22);S24=S1(24);

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

W11=10*X3*r1^6*L6;             W12=7*sqrt(15)*X3*r1^8*L8;      W13=12*sqrt(21)*X3*r1^10*L10;    W14=55*sqrt(27)*X3*r1^12*L12;     W15=26*sqrt(33)*X3*r1^14*L14;
W21=7*sqrt(15)*X5*r1^8*L8;     W22=126*X5*r1^10*L10;           W23=66*sqrt(35)*X5*r1^12*L12;    W24=429*sqrt(5)*X5*r1^14*L14;     W25=273*sqrt(55)*X5*r1^16*L16;
W31=12*sqrt(21)*X7*r1^10*L10;  W32=66*sqrt(35)*X7*r1^12*L12;   W33=1716*X7*r1^14*L14;           W34=2145*sqrt(7)*X7*r1^16*L16;    W35=1768*sqrt(77)*X7*r1^18*L18; 
W41=55*sqrt(27)*X9*r1^12*L12;  W42=429*sqrt(5)*X9*r1^14*L14;   W43=2145*sqrt(7)*X9*r1^16*L16;   W44=24310*X9*r1^18*L18;           W45=25194*sqrt(11)*X9*r1^20*L20; 
W51=26*sqrt(33)*X11*r1^14*L14; W52=273*sqrt(55)*X11*r1^16*L16; W53=1768*sqrt(77)*X11*r1^18*L18; W54=25194*sqrt(11)*X11*r1^20*L20; W55=352716*X11*r1^22*L22;

W=[W11, W12, W13, W14, W15; W21, W22, W23, W24, W25;W31, W32, W33, W34, W35;...
       W41, W42, W43, W44, W45; W51, W52, W53, W54, W55];

AP2=eye(10)+W;
PP=-X1*N1*inv(AP2)*N2;  
JK=J+PP; % Matriz Z
Z=JK;
RR5 =coeficientes_efectivos(cm,x1,x2,V3,V2,Z);

end


% %**************************************************************************
% % FORMULA DE LOS COEFICIENTES EFECTIVOS (PROBLEMAS Antiplano)
% %**************************************************************************
function RR =coeficientes_efectivos(cm,x1,x2,V3,V2,Z)

Delta= (V2*(1-x1)*(x1+x2)+2*V3*x1*(1-x2)) / ...
         ((1+x1)*(x1+x2)*(V2+V3)+(1-x1)*(x1-x2)*(V3));
 %F=(1-x1+2*V3*(x1-x2)*x1/(V2*(x1+x2)+2*V3*x1));
% G=1/(V2*(x1+x2)+2*V3*x1);
DtZ=det(Z);
c55=cm*(1-2*(V2+V3)*Delta*(Z(2,2))/(DtZ));
  c45=2*cm*((V2+V3)*Delta*(Z(2,1))/(DtZ));
  c54=2*((V2+V3)*Delta*(Z(1,2))/(DtZ));
c44=cm*(1-2*(V2+V3)*Delta*(Z(1,1))/(DtZ));

      
RR=[c55,c44,c54,c45];
end

