clearvars
clear all

% tabla Table 4 Jiang 2004
%Effective elastic moduli C44=G, C55=G, C45=G (G'/G 120)
Jiang2004 =[1.21780, 0.00133, 1.21852, 1.22306, -0.00471, 1.21365;
           1.48816, 0.00672, 1.49176, 1.51608, -0.02383 ,1.46841;
           1.83374, 0.01968, 1.84428, 1.92172, -0.07065, 1.78042;
           2.29477, 0.04739, 2.32016, 2.53343, -0.17671, 2.18002;
           2.95203, 0.10602, 3.00884, 3.62069, -0.44144, 2.73782;
           4.00126, 0.23812, 4.12884, 0,0,0;
           0,0,0,0,0,0;
           0,0,0,0,0,0;
           0,0,0,0,0,0];

 %%%%%%%%%%%%%%%%%%%%%
% 
load S10_45_75_90_delta1_2; 

cm=1;
cf=120;
x=cf/cm;
X=(1-x)/(1+x);
xas=x;

tt=0;
for ii=[5,9]
     tt=tt+1;
% ************************************************************************
S1=S10_45_75_90_delta1_2(:,ii);
 delta1=2*ZZ(ii);
 Periodos2=(cos(theta(ii))+1i*sin(theta(ii)));
 w1=1;
 w2=Periodos2;
rm=min(sin(theta(ii)/2),0.5);
percolacion(tt) =pi*rm.^2./sin(theta(ii));

  
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

if ii==5
lam=[0:.0125:.65]; 
Lon(1)=length(lam);
else
lam=[0:.0125:.8, .81]; 
Lon(2)=length(lam);
end
% *************************************************************************

for j=1:length(lam)
  V2=lam(j);
  r = sqrt(V*V2/pi);  % *** Radio de la fibra más la mesophase ************************************

 
 J=eye(2)+X*r^2*[(H1+h12),(-h22);(-h22),(H1-h12)];
 
 
 CENo0(j,:,tt)=coeficientes_efectivos2f(cm,V2,X,J); % VanFo Fi
 CENo1(j,:,tt) =ahm2f_Paralelogramo_Cortas_2023_No1(cm,X,V2,S1,r,J);
 CENo2(j,:,tt) =ahm2f_Paralelogramo_Cortas_2023_No2(cm,X,V2,S1,r,J); 
 CENo3(j,:,tt) =ahm2f_Paralelogramo_Cortas_2023_No3(cm,X,V2,S1,r,J); 
 CENo4(j,:,tt) =ahm2f_Paralelogramo_Cortas_2023_No4(cm,X,V2,S1,r,J); 
 CENo5(j,:,tt) =ahm2f_Paralelogramo_Cortas_2023_No5(cm,X,V2,S1,r,J); 
 CENoGm(j,:,tt) =conductividad_efectiva_orden_m(cm,X,V2,S1,r,J,30);
                                 
end
    
angulo(tt,:)=['$\theta=$', num2str(theta(ii)*180/pi,3),'$^o$'];
end


for tt=1:2
matriz=[ '$\kappa_{1}=$',num2str(cm,5)];
fibra=[ '$\kappa_{2}=$',num2str(cf,5)];
Modulo=[ '|w2|= ',num2str(abs(w2),5)];
roo=[ '$\rho=$',num2str(cf/cm,5)];

if tt==1
lam=[0:.0125:.65]; 
Lon(1)=length(lam);
t2=[9:8:Lon(1)-4,Lon(1)];
else
lam=[0:.0125:.8, .81]; 
Lon(2)=length(lam);
t2=[9:8:Lon(2)-1,Lon(2)];
end

%textoY=['$\hat{\kappa}_{11}/\kappa_1$';'$\hat{\kappa}_{22}/\kappa_1$';'$\hat{\kappa}_{12}/\kappa_1$'];
disp(Modulo)
disp(angulo(tt,:))
disp(matriz)
disp(fibra)
Tabla=zeros(Lon(tt),7);
Tabla(:,:)=[CENo0(1:Lon(tt),1,tt), CENo1(1:Lon(tt),1,tt),CENo2(1:Lon(tt),1,tt),CENo3(1:Lon(tt),1,tt),...
            CENo4(1:Lon(tt),1,tt),CENo5(1:Lon(tt),1,tt), CENoGm(1:Lon(tt),1,tt)];


 disp('Vf   	\hat{\kappa_11}/\kappa_1 ') 
 disp('V2          o_0,         o_1,        o_2,         o_3,         o_4,        o_5           o_30         Jiang')
 disp(' -------------------------------------------')
 disp(num2str([lam(t2)', Tabla(t2,:) Jiang2004(1:length(t2),6-3*(tt-1))],6))
 
 disp(' -------------------------------------------')
 Tabla(:,:)=[CENo0(1:Lon(tt),2,tt), CENo1(1:Lon(tt),2,tt),CENo2(1:Lon(tt),2,tt),CENo3(1:Lon(tt),2,tt),...
            CENo4(1:Lon(tt),2,tt),CENo5(1:Lon(tt),2,tt),CENoGm(1:Lon(tt),2,tt)];
  disp('Vf   	\hat{\kappa_22}/\kappa_1 ')  
  disp('V2          o_0,         o_1,        o_2,         o_3,         o_4,        o_5           o_30         Jiang')
 disp(' -------------------------------------------')
disp(num2str([lam(t2)', Tabla(t2,:), Jiang2004(1:length(t2),4-3*(tt-1))],6))
 
 Tabla(:,:)=[CENo0(1:Lon(tt),3,tt), CENo1(1:Lon(tt),3,tt),CENo2(1:Lon(tt),3,tt),CENo3(1:Lon(tt),3,tt),...
            CENo4(1:Lon(tt),3,tt),CENo5(1:Lon(tt),3,tt),CENoGm(1:Lon(tt),3,tt)];
 disp(' -------------------------------------------')
  disp('Vf   	\hat{\kappa_12}/\kappa_1 ')    
 disp('V2          o_0,         o_1,        o_2,         o_3,         o_4,        o_5           o_30         Jiang')
 disp(' -------------------------------------------')
disp(num2str([lam(t2)', Tabla(t2,:), Jiang2004(1:length(t2),5-3*(tt-1))],5))
end

% colores=[[0,0,1];[.6,.6,0];[0,1,0];[1,0,0];[0,1,1];...
%     [1,0,1];[.8,0.5,.8];[0.5,0.5,0.5];[0.5,0.5,1];[0.8,0.8,0];[0.1,0.8,0];[.8,0,0];[0,.8,.8]];
% marca=['ko-'; 'ko-';'ko-';'ko-'; 'ko-';'ko-';'ko-';'ks-'];
% letras=['a)';'b)';'c)';'d)';'e)';'f)']  ;
% figure(1)
%   
% for tt=1:2
%   
% for jj=1:3
%  subplot(2,3,3*(tt-1)+jj)
% hold on
% plot(lam(1:Lon(tt)),CENoGm(1:Lon(tt),jj,tt) ,marca(8,:),'LineWidth',1,...
%                        'MarkerEdgeColor',colores(8,:),...
%                        'MarkerFaceColor','w',...
%                         'MarkerSize',6)
% plot(lam(1:Lon(tt)),CENo5(1:Lon(tt),jj,tt) ,marca(6,:),'LineWidth',1,...
%                        'MarkerEdgeColor',colores(6,:),...
%                        'MarkerFaceColor',colores(6,:),...
%                         'MarkerSize',3)
% 
% 
% plot(lam(1:Lon(tt)),CENo4(1:Lon(tt),jj,tt) ,marca(5,:),'LineWidth',1,...
%                        'MarkerEdgeColor',colores(5,:),...
%                        'MarkerFaceColor',colores(5,:),...
%                         'MarkerSize',3)
% 
% plot(lam(1:Lon(tt)),CENo3(1:Lon(tt),jj,tt) ,marca(4,:),'LineWidth',1,...
%                        'MarkerEdgeColor',colores(4,:),...
%                        'MarkerFaceColor',colores(4,:),...
%                         'MarkerSize',3) 
% plot(lam(1:Lon(tt)),CENo2(1:Lon(tt),jj,tt) ,marca(3,:),'LineWidth',1,...
%                        'MarkerEdgeColor',colores(3,:),...
%                        'MarkerFaceColor',colores(3,:),...
%                         'MarkerSize',3) 
% plot(lam(1:Lon(tt)),CENo1(1:Lon(tt),jj,tt) ,marca(2,:),'LineWidth',1,...
%                        'MarkerEdgeColor',colores(2,:),...
%                        'MarkerFaceColor',colores(2,:),...
%                         'MarkerSize',3) 
% plot(lam(1:Lon(tt)),CENo0(1:Lon(tt),jj,tt) ,marca(1,:),'LineWidth',1,...
%                        'MarkerEdgeColor',colores(1,:),...
%                        'MarkerFaceColor',colores(1,:),...
%                         'MarkerSize',3)
% 
% hold off
% 
% if tt==2
% xlabel('Fiber volume fraction','fontsize',12,'interpreter','latex')
% end
% xlim([0, lam(Lon(tt))+.1])
% xticks([0:.1:lam(Lon(tt))+.1]);
% text(.2,CENoGm(Lon(tt),jj,tt),[letras(3*(tt-1)+jj,:)],'fontsize',14,'interpreter','latex')
% set(gca, 'FontSize', 12)
% ylabel(textoY(jj,:),'fontsize',18,'interpreter','latex')
% grid on
% grid minor
% box on
% % title([roo, angulo(tt,:)],'fontsize',18,'interpreter','latex')
% title( [angulo(tt,:)],'fontsize',12,'interpreter','latex')
%  if jj==1 & tt==1
% % ylim([0,floor(max(CENoG(:,1,tt))+1)])
% legend('$O_{30}$','$O_5$','$O_4$','$O_3$','$O_2$','$O_1$','$O_0$','fontsize',14,'Location','northwest','interpreter','latex')
% text(.1,2.5,['$\rho=$',num2str(xas,3)],'fontsize',14,'interpreter','latex')
% % text(.01,9,fibra, 'fontsize',14,'interpreter','latex') 
% % text(.01,8,angulo(tt,:), 'fontsize',12,'interpreter','latex')
%  end
% end
% end


% Funciones auxiliares
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

