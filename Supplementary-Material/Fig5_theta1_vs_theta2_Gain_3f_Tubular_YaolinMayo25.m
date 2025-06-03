%Este programa resuelve el problema elástico trifásico 
%       DIRECTAMENTE.

%INTRODUCCION DE LOS DATOS DE LA MATRIZ   

clearvars;close all;
clear all
% %%%%% ejes no principales
load S10_45_75_90_delta1_2;

cm=1;
cf=cm;
cmf=500;
tt=0;

ii=[5,9];
V2masV3a=[.65 ,.6,.5;.784 ,.7,.65];
for tt=1:2
 V2masV3=[];    
Sk=S10_45_75_90_delta1_2(:,ii(tt));
Periodos=cos(theta(ii(tt)))+1i*sin(theta(ii(tt))); % w2
w2=Periodos;
w1=1;

delta1=2*ZZ(ii(tt));
delta2=2*ZZ1(ii(tt));   %delta1Leg.*Periodos-2*pi*1i;  

V = abs(w1)*abs(w2)*sin(theta(ii(tt))); % *** Volumen de la celda periodica *******
% *** expresiones de H1 y H2 en el problema I23 ********
H1=(conj(delta1)*conj(w2)-conj(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
H2=((delta1)*conj(w2)-(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
% *************************************************************************

rm=min(sin(theta(ii(tt))/2),0.5);
percolacion =pi*rm.^2./sin(theta(ii(tt)));

%% Selecc volumen fiber fraction
if tt==1
    V2masV3=V2masV3a(1,:);
else
V2masV3=V2masV3a(2,:);
end
% V2masV3=[.5];
%% Determinar fraccion de volumen  de la mesofase y la fibbra central
for jj=1:length(V2masV3)
 alpha=[1e-5, .001:.005:.1, .125:.025:.95,.96,.97,.98,.99, 1];   
me(jj,:)=alpha.*V2masV3(jj);
% lam2(jj,:)=;
lam2(jj,:)=V2masV3(jj)-me(jj,:);
end

%% calcular las propiedades efectivas de tres fases
for kk=1:length(V2masV3)

 for jj=1:length(lam2(kk,:))
        
 RR1= ahm3f_Paralelogramo_Cortas_2023_No4(cm,cf,cmf,lam2(kk,jj),me(kk,jj),V,Sk,H1,H2);
    
 RR=[(RR1(1)),(RR1(2)),(RR1(3)),  (RR1(4)) ];
 
 CETAHM_3f(jj,:,kk,tt)=[RR]; 
 end
 
end
%%% calcular las propiedades efectivas de dos fases
 
for kk=1:length(V2masV3)
 for jj=1:length(me(kk,:))
        
  RR1=ahm2f_Paralelogramo_Cortas_2023_No4(cm,cmf,me(kk,jj),V,Sk,H1,H2);
     
  RR=[(RR1(1)),(RR1(2)),(RR1(3)),  (RR1(4)) ];
 
 CETAHM_2f(jj,:,kk,tt)=[RR]; % Para obtener la matriz de las posisiones del tensor por vf y ángulo
  end
 
end

end

  elasticos=['$\hat \kappa_{11}$';'$\hat \kappa_{22}$'; '$\hat \kappa_{12}$'; '$\hat \kappa_{21}$' ];
  
 matriz=[ 'p_1= ',num2str(cm,5)];
fibra=[ 'p_2= ',num2str(cf,5)];
mesofase=[ 'p_I= ',num2str(cmf,5)];

% angulo(:,tt)=['angle of the periodic cell \theta= ' num2str(theta(ii(tt))*180/pi,3),'^o']; 
  
 disp(['matrix =',num2str(cm)])
  disp(['fiber=',num2str(cm)])
disp(['Ring =',num2str(cmf)])

disp(['angle of the periodic cell= ' num2str(theta(ii(tt))*180/pi,3)])
% 
%     disp('******************************************')
%     disp([' h ,   Real(kappa)AHM , Imag(kappa)AHM,  ' ])
% disp( num2str([h', CETAHM_3f(:,1), CETAHM_3f(:,2), CETAHM_3f(:,3) CETAHM_3f(:,4), CETAHM_3f(:,5), CETAHM_3f(:,6)],5))
letras=['a)';'b)';'c)';'d)';'e)';'f)']  ;
colores=[[0,0,1];[.6,.6,0];[0,1,0];[1,0,0];[0,1,1];[1,0,1];[0,0,0];[0.5,0.5,0.5];[0.5,0.5,1];[0.8,0.8,0];[0.1,0.8,0];[.8,0,0];[0,.8,.8]];
Gain=['$\kappa^{gain}_{11}$';'$\kappa^{gain}_{22}$';'$\kappa^{gain}_{12}$'];
%%%%%%%%%%%%%%%%%%%%%%%

figure(1);
 for tt=1:2
   for jj=1:3
      subplot(2,3,3*(tt-1)+jj)
     for kk=1:length(V2masV3a(tt,:))
     hold on
     % Seleccionar estilo basado en el índice kk
      switch mod(kk,8)
        case 1
            linestyle = '-';  linewidth = 1.5;
        case 2
            linestyle = '--'; linewidth = 1.5;
        case 3
            linestyle = ':';  linewidth = 2.0; 
        case 4
            linestyle = '-';  linewidth = 2.0;
       end
       h(kk) =plot(alpha,CETAHM_3f(:,jj,kk,tt)./CETAHM_2f(:,jj,kk,tt),...
                      'LineStyle', linestyle, ...
                       'LineWidth', linewidth, ...
                       'Color', colores(kk,:));

      legendEntries{kk} = ['\textbf{$V_2+V_3=$}',num2str(V2masV3a(tt,kk),3)];  
     end
         h(4) =plot(alpha,CETAHM_2f(:,jj,kk,tt)./CETAHM_2f(:,jj,kk,tt),...
                      'LineStyle', linestyle, ...
                       'LineWidth', linewidth, ...
                       'Color', colores(4,:));
         legendEntries{4} = ['\textbf{$V_3=$}',num2str(0)];  
       hold off
     
ylabel( Gain(jj,:),'fontsize',18,'interpreter','latex')
ha=max(CETAHM_3f(:,jj,1,tt)./CETAHM_2f(:,jj,1,tt))+1;
ylim([0,ha ])
text(.1,ha-1,[letras(3*(tt-1)+jj,:)],'fontsize',14,'interpreter','latex')
title(['\textbf{$\theta=$}',num2str(theta(ii(tt))*180/pi,3),'$^o$',...
       ', ','\textbf{$\rho=$}',num2str(cmf/cm,3)],'fontsize',12,'interpreter','latex')

    if jj==3 
    % Leyenda
    legend(h, legendEntries, 'FontSize', 12, 'Location', 'best', 'interpreter', 'latex', 'FontWeight', 'bold','NumColumns',7);
    end
   if tt==2
       xlabel('$\alpha=\frac{V_2}{V_2+V_3}$','fontsize',14,'interpreter','latex')
   end
grid on
grid minor
box on
   end
 end
 
 %%% orden del sistema No=4
function RR3 =ahm3f_Paralelogramo_Cortas_2023_No4(cm,cf,cmf,V3,V2,V,S1,H1,H2)
x1=cmf/cm;
x2=cf/cm;
r1 = sqrt(V*(V2+V3)/pi) ;%% radio externo
X1=((1-x1)*(x1+x2)*(V2+V3)+(1+x1)*(x1-x2)*(V3))/...
             ((1+x1)*(x1+x2)*(V2+V3)+(1-x1)*(x1-x2)*(V3));
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
J=eye(2)+X1*r1^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)];
              

S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);S16=S1(16);S18=S1(18);

X3=((1-x1)*(x1+x2)*(V2+V3)^3+(1+x1)*(x1-x2)*(V3)^3)/...
             ((1+x1)*(x1+x2)*(V2+V3)^3+(1-x1)*(x1-x2)*(V3)^3);   
         
X5=((1-x1)*(x1+x2)*(V2+V3)^5+(1+x1)*(x1-x2)*(V3)^5)/...
             ((1+x1)*(x1+x2)*(V2+V3)^5+(1-x1)*(x1-x2)*(V3)^5);
X7=((1-x1)*(x1+x2)*(V2+V3)^7+(1+x1)*(x1-x2)*(V3)^7)/...
             ((1+x1)*(x1+x2)*(V2+V3)^7+(1-x1)*(x1-x2)*(V3)^7);
X9=((1-x1)*(x1+x2)*(V2+V3)^9+(1+x1)*(x1-x2)*(V3)^9)/...
             ((1+x1)*(x1+x2)*(V2+V3)^9+(1-x1)*(x1-x2)*(V3)^9);


N1=[sqrt(3)*r1^4*real(S4),-sqrt(3)*r1^4*imag(S4),sqrt(5)*r1^6*real(S6),-sqrt(5)*r1^6*imag(S6), sqrt(7)*r1^8*real(S8),-sqrt(7)*r1^8*imag(S8),3*r1^10*real(S10),-3*r1^10*imag(S10);
    -sqrt(3)*r1^4*imag(S4), -sqrt(3)*r1^4*real(S4),-sqrt(5)*r1^6*imag(S6), -sqrt(5)*r1^6*real(S6),-sqrt(7)*r1^8*imag(S8), -sqrt(7)*r1^8*real(S8),-3*r1^10*imag(S10), -3*r1^10*real(S10)];
N2=[X3*sqrt(3)*r1^4*real(S4),-X3*sqrt(3)*r1^4*imag(S4),  X5*sqrt(5)*r1^6*real(S6),-X5*sqrt(5)*r1^6*imag(S6),X7*sqrt(7)*r1^8*real(S8),-X7*sqrt(7)*r1^8*imag(S8),3*X9*r1^10*real(S10),-3*X9*r1^10*imag(S10);
    -X3*sqrt(3)*r1^4*imag(S4), -X3*sqrt(3)*r1^4*real(S4),-X5*sqrt(5)*r1^6*imag(S6), -X5*sqrt(5)*r1^6*real(S6),-X7*sqrt(7)*r1^8*imag(S8), -X7*sqrt(7)*r1^8*real(S8),-3*X9*r1^10*imag(S10), -3*X9*r1^10*real(S10)];
N2=N2';

W=[10*r1^6*X3*real(S6),-10*r1^6*X3*imag(S6),7*sqrt(15)*r1^8*X3*real(S8),-7*sqrt(15)*r1^8*X3*imag(S8),12*sqrt(21)*r1^10*X3*real(S10),-12*sqrt(21)*r1^10*X3*imag(S10),55*sqrt(27)*r1^12*X3*real(S12),-55*sqrt(27)*r1^12*X3*imag(S12);...
    -10*r1^6*X3*imag(S6),-10*r1^6*X3*real(S6),-7*sqrt(15)*r1^8*X3*imag(S8),-7*sqrt(15)*r1^8*X3*real(S8),-12*sqrt(21)*r1^10*X3*imag(S10),-12*sqrt(21)*r1^10*X3*real(S10),-55*sqrt(27)*r1^12*X3*imag(S12),-55*sqrt(27)*r1^12*X3*real(S12);...
    7*sqrt(15)*r1^8*X5*real(S8),-7*sqrt(15)*r1^8*X5*imag(S8),126*r1^10*X5*real(S10),-126*r1^10*X5*imag(S10),66*sqrt(35)*r1^12*X5*real(S12),-66*sqrt(35)*r1^12*X5*imag(S12),429*sqrt(5)*r1^14*X5*real(S14),-429*sqrt(5)*r1^14*X5*imag(S14);...
    -7*sqrt(15)*r1^8*X5*imag(S8),-7*sqrt(15)*r1^8*X5*real(S8),-126*r1^10*X5*imag(S10),-126*r1^10*X5*real(S10),-66*sqrt(35)*r1^12*X5*imag(S12),-66*sqrt(35)*r1^12*X5*real(S12),-429*sqrt(5)*r1^14*X5*imag(S14),-429*sqrt(5)*r1^14*X5*real(S14);...
    12*sqrt(21)*r1^10*X7*real(S10),-12*sqrt(21)*r1^10*X7*imag(S10),66*sqrt(35)*r1^12*X7*real(S12),-66*sqrt(35)*r1^12*X7*imag(S12),1716*r1^14*X7*real(S14),-1716*r1^14*X7*imag(S14),2145*sqrt(7)*r1^16*X7*real(S16),-2145*sqrt(7)*r1^16*X7*imag(S16);...
    -12*sqrt(21)*r1^10*X7*imag(S10),-12*sqrt(21)*r1^10*X7*real(S10),-66*sqrt(35)*r1^12*X7*imag(S12),-66*sqrt(35)*r1^12*X7*real(S12),-1716*r1^14*X7*imag(S14),-1716*r1^14*X7*real(S14),-2145*sqrt(7)*r1^16*X7*imag(S16),-2145*sqrt(7)*r1^16*X7*real(S16);...
      55*sqrt(27)*r1^12*X9*real(S12), -55*sqrt(27)*r1^12*X9*imag(S12),429*sqrt(5)*r1^14*X9*real(S14),-429*sqrt(5)*r1^14*X9*imag(S14), 2145*sqrt(7)*r1^16*X9*real(S16),-2145*sqrt(7)*r1^16*X9*imag(S16),24310*r1^18*X9*real(S18),-24310*r1^18*X9*imag(S18);...
      -55*sqrt(27)*r1^12*X9*imag(S12), -55*sqrt(27)*r1^12*X9*real(S12),-429*sqrt(5)*r1^14*X9*imag(S14),-429*sqrt(5)*r1^14*X9*real(S14),-2145*sqrt(7)*r1^16*X9*imag(S16),-2145*sqrt(7)*r1^16*X9*real(S16),-24310*r1^18*X9*imag(S18),-24310*r1^18*X9*real(S18)];
     
AP2=eye(8)+W;
PP=-X1*N1*inv(AP2)*N2;  
JK1=J+PP; % Matriz Z
Z=JK1;

RR3 =coeficientes_efectivos3f(cm,x1,x2,V3,V2,X1,Z);

end

%%% orden del sistema No=4
function RR3 =ahm2f_Paralelogramo_Cortas_2023_No4(cm,cf,V2,V,S1,H1,H2)
 
x=cf/cm;
r1 = sqrt(V*(V2)/pi) ;%% radio externo
X=(1-x)/(1+x);
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
J=eye(2)+X*r1^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)];

S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);S16=S1(16);S18=S1(18);

N1=[sqrt(3)*r1^4*real(S4),-sqrt(3)*r1^4*imag(S4),sqrt(5)*r1^6*real(S6),-sqrt(5)*r1^6*imag(S6), sqrt(7)*r1^8*real(S8),-sqrt(7)*r1^8*imag(S8),3*r1^10*real(S10),-3*r1^10*imag(S10);
    -sqrt(3)*r1^4*imag(S4), -sqrt(3)*r1^4*real(S4),-sqrt(5)*r1^6*imag(S6), -sqrt(5)*r1^6*real(S6),-sqrt(7)*r1^8*imag(S8), -sqrt(7)*r1^8*real(S8),-3*r1^10*imag(S10), -3*r1^10*real(S10)];
N2=N1';

W=[10*r1^6*real(S6),-10*r1^6*imag(S6),7*sqrt(15)*r1^8*real(S8),-7*sqrt(15)*r1^8*imag(S8),12*sqrt(21)*r1^10*real(S10),-12*sqrt(21)*r1^10*imag(S10),55*sqrt(27)*r1^12*real(S12),-55*sqrt(27)*r1^12*imag(S12);...
    -10*r1^6*imag(S6),-10*r1^6*real(S6),-7*sqrt(15)*r1^8*imag(S8),-7*sqrt(15)*r1^8*real(S8),-12*sqrt(21)*r1^10*imag(S10),-12*sqrt(21)*r1^10*real(S10),-55*sqrt(27)*r1^12*imag(S12),-55*sqrt(27)*r1^12*real(S12);...
    7*sqrt(15)*r1^8*real(S8),-7*sqrt(15)*r1^8*imag(S8),126*r1^10*real(S10),-126*r1^10*imag(S10),66*sqrt(35)*r1^12*real(S12),-66*sqrt(35)*r1^12*imag(S12),429*sqrt(5)*r1^14*real(S14),-429*sqrt(5)*r1^14*imag(S14);...
    -7*sqrt(15)*r1^8*imag(S8),-7*sqrt(15)*r1^8*real(S8),-126*r1^10*imag(S10),-126*r1^10*real(S10),-66*sqrt(35)*r1^12*imag(S12),-66*sqrt(35)*r1^12*real(S12),-429*sqrt(5)*r1^14*imag(S14),-429*sqrt(5)*r1^14*real(S14);...
    12*sqrt(21)*r1^10*real(S10),-12*sqrt(21)*r1^10*imag(S10),66*sqrt(35)*r1^12*real(S12),-66*sqrt(35)*r1^12*imag(S12),1716*r1^14*real(S14),-1716*r1^14*imag(S14),2145*sqrt(7)*r1^16*real(S16),-2145*sqrt(7)*r1^16*imag(S16);...
    -12*sqrt(21)*r1^10*imag(S10),-12*sqrt(21)*r1^10*real(S10),-66*sqrt(35)*r1^12*imag(S12),-66*sqrt(35)*r1^12*real(S12),-1716*r1^14*imag(S14),-1716*r1^14*real(S14),-2145*sqrt(7)*r1^16*imag(S16),-2145*sqrt(7)*r1^16*real(S16);...
      55*sqrt(27)*r1^12*real(S12), -55*sqrt(27)*r1^12*imag(S12),429*sqrt(5)*r1^14*real(S14),-429*sqrt(5)*r1^14*imag(S14), 2145*sqrt(7)*r1^16*real(S16),-2145*sqrt(7)*r1^16*imag(S16),24310*r1^18*real(S18),-24310*r1^18*imag(S18);...
      -55*sqrt(27)*r1^12*imag(S12), -55*sqrt(27)*r1^12*real(S12),-429*sqrt(5)*r1^14*imag(S14),-429*sqrt(5)*r1^14*real(S14),-2145*sqrt(7)*r1^16*imag(S16),-2145*sqrt(7)*r1^16*real(S16),-24310*r1^18*imag(S18),-24310*r1^18*real(S18)];
     
AP2=eye(8)+X*W;
PP=-X^2*N1*inv(AP2)*N2;  
JK=J+PP; % Matriz Z
Z=JK;

RR3 =coeficientes_efectivos2f(cm,V2,X,Z);

end



function RRG = ahm3f_elastico_p13_p23_Paralelogramo_Agosto(cm,cf,cmf,V3,V2,V,S1,H1,H2,No)
% *************************************************************************
% *** Metodo de Homogenizacion Asintotica *********************************
% *************************************************************************
% ***   ***   ***   Simetría rómbica. Problema ANTIPLANO   ***   ***   ****
% ***   ***   ***   Constituyentes Piezoelectricas (TIM)   ***   ***   ****
% ***   ***   ***   ***   ***   Conatcto Imperfecto   **   ***   ***   **** 
% *************************************************************************
S=S1;
%No=5; % *** Orden del sistema ******************************************

x1=cmf/cm;
x2=cf/cm;

% % *** expresiones de H1 y H2 en el problema I23 ********
% H1=(conj(delta1)*conj(w2)-conj(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
% H2=((delta1)*conj(w2)-(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
% % *************************************************************************
 
r1 = sqrt(V*(V2+V3)/pi);  % *** Radio de la fibra más la mesophase ************************************

%Matriz Wnk
for n=1:2:2*No+1,
   for k=1:2:35,

Wnk(n,k)=(S(n+k))*(r1^(n+k))*((factorial(n+k-1)/(factorial(n-1)*factorial(k-1)))/...
               sqrt(n*k));

   end
end
RWnk=real(Wnk);
IWnk=imag(Wnk);
di=eye(20); 
% *************************************************************************
for s=1:2:2*No+1,

       Xp(s)=((1-x1)*(x1+x2)*(V2+V3)^s+(1+x1)*(x1-x2)*(V3)^s)/...
             ((1+x1)*(x1+x2)*(V2+V3)^s+(1-x1)*(x1-x2)*(V3)^s);
       
end


%      %*************** MATRIZ del Sistema Problema 13 y 23 . Problema antiplano
  for t=1:No,
   for s=1:No ;

   AP(:,:,t,s)=di(t,s)*[1,0;0,1] + Xp(2*t+1)*[RWnk(2*s+1,2*t+1), -IWnk(2*s+1,2*t+1);...
                                        -IWnk(2*s+1,2*t+1), -RWnk(2*s+1,2*t+1)];
  
   end
  end
  
  for t=1:No
      
    NP(:,:,t)=Xp(2*t+1)*[RWnk(1,2*t+1), -IWnk(1,2*t+1);...
                         -IWnk(1,2*t+1),-RWnk(1,2*t+1)];
                     
    NPa(:,:,t)=[RWnk(2*t+1,1), -IWnk(2*t+1,1);...
                         -IWnk(2*t+1,1),-RWnk(2*t+1,1)];
     
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
 
  
  PP=-Xp(1)*NP1a*inv(AP2)*NP1; 
  

h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);

J=[1+Xp(1)*r1^2*(h11+h12),Xp(1)*r1^2*(h21-h22);...
     -Xp(1)*r1^2*(h21+h22),1+Xp(1)*r1^2*(h11-h12)]; % cortas
 
 
JK=J+PP;
Z=JK
DtZ=det(Z);
% E=(V2*(x1+x2)*(1-x1)+2*V3*(x1)*(1-x2))/...
%     (V2*(x1+x2)*(1+x1)+2*V3*(x1)*(1+x2));
E=Xp(1);
%**************************************************************************
% FORMULA DE LOS COEFICIENTES EFECTIVOS (PROBLEMAS Antiplano)
%**************************************************************************
RRG =coeficientes_efectivos3f(cm,x1,x2,V3,V2,Xp(1),Z);

end

function RR =coeficientes_efectivos2f(cm,V2,X,Z) 

DtZ=det(Z);
c55=cm*(1-2*V2*X*Z(2,2)/DtZ);
c45=cm*(2*V2*X*Z(2,1)/DtZ);
c54=cm*(2*V2*X*Z(1,2)/DtZ);
c44=cm*(1-2*V2*X*Z(1,1)/DtZ);
        
RR=[c55,c44,c54,c45];
end
%**************************************************************************
% FORMULA DE LOS COEFICIENTES EFECTIVOS 3Fases (PROBLEMAS Antiplano)
%**************************************************************************
function RR =coeficientes_efectivos3f(cm,x1,x2,V3,V2,X1,Z)

Delta=(V2*(1-x1)*(x1+x2)+2*V3*x1*(1-x2)) / ...
         ((1+x1)*(x1+x2)*(V2+V3)+(1-x1)*(x1-x2)*(V3));
DtZ=det(Z);
c55=cm*(1-2*(V2+V3)*Delta*(Z(2,2))/(DtZ));
  c45=2*cm*((V2+V3)*Delta*(Z(2,1))/(DtZ));
  c54=2*cm*((V2+V3)*Delta*(Z(1,2))/(DtZ));
c44=cm*(1-2*(V2+V3)*Delta*(Z(1,1))/(DtZ));

      
RR=[c55,c44,c54,c45];
end


