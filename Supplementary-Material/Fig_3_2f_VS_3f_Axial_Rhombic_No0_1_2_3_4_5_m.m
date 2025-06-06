clearvars

ii=5; % %% 45
lam=[.01, .025:.025:.60];
% lam=[0, .1:.1:.6 ];  
Lon=length(lam);

 %%%%%%%%%%%%%%%%%%%%%
load S10_45_75_90_delta1_2; 

 S1=S10_45_75_90_delta1_2(:,ii);
 delta1=2*ZZ(ii);
 Periodos2=(cos(theta(ii))+1i*sin(theta(ii)));
 w1=1;
 w2=Periodos2;

delta1=delta1;
delta2=delta1.*w2-2*pi*1i; 

V = abs(w1)*abs(w2)*sin(theta(ii)); % *** Volumen de la celda periodica *******

H1=(conj(delta1)*conj(w2)-conj(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
H2=((delta1)*conj(w2)-(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
cm=1;
cmf=12;
cf=120;
x1=cmf/cm;
x2=cf/cm;
% *************************************************************************
me=[0,5e-2]; %% interfase volume fraction
xas=[x1,x2];

    
for tt=1:2
for j=1:length(lam)
    
    V3=lam(j);
    V2=me(tt);
  r1 = sqrt(V*(V2+V3)/pi);  % *** Radio de la fibra m�s la mesophase ************************************
X1=((1-x1)*(x1+x2)*(V2+V3)+(1+x1)*(x1-x2)*(V3))/...
             ((1+x1)*(x1+x2)*(V2+V3)+(1-x1)*(x1-x2)*(V3));
 J=eye(2)+X1*r1^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)];
 
 CENo0(j,:,tt)=coeficientes_efectivos(x1,x2,V3,V2,J); 
 CENo1(j,:,tt) =ahm3f_Paralelogramo_Cortas_2025_No1(x1,x2,V3,V2,S1,r1,J);
 CENo2(j,:,tt) =ahm3f_Paralelogramo_Cortas_2025_No2(x1,x2,V3,V2,S1,r1,J); 
 CENo3(j,:,tt) =ahm3f_Paralelogramo_Cortas_2025_No3(x1,x2,V3,V2,S1,r1,J);
 CENo4(j,:,tt) =ahm3f_Paralelogramo_Cortas_2025_No4(x1,x2,V3,V2,S1,r1,J);
 CENo5(j,:,tt) =ahm3f_Paralelogramo_Cortas_2025_No5(x1,x2,V3,V2,S1,r1,J);
 CENo6(j,:,tt) =ahm3f_Paralelogramo_Cortas_2025_No6(x1,x2,V3,V2,S1,r1,J);
 CENoG(j,:,tt) = conductividad_efectiva_orden_m(cm,cf,cmf,V3,V2,V,S1,r1,J,6);
                 
end
end


matriz=[ 'cm= ',num2str(cm,5)];
fibra=[ 'cf= ',num2str(cf,5)];
mesofase=[ 'cmf= ',num2str(cmf,5)];
rho1=[ '$\rho_1$','=',num2str(cmf/cm,5),];
rho2=[ '$\rho_2$','=',num2str(cf/cm,5)];
angulo=['$theta=$ ' num2str(theta(ii)*180/pi,3),'$^o$'];
Modulo=[ '|w2|= ',num2str(abs(w2),5)];
volmesofase=[ '$V_2= $',num2str(V2,3)];
disp(Modulo)
disp(angulo)
disp(matriz)
disp(fibra)
disp(mesofase)

textoY=['$\hat{\kappa}_{11}/\kappa_1$';'$\hat{\kappa}_{22}/\kappa_1$';'$\hat{\kappa}_{12}/\kappa_1$'];
%%% Dos  fases

for tt=1:2

Vmf=['$V_{2}=$',num2str(num2str(me(tt)),3)];

matriz=[ 'cm= ',num2str(cm,5)];
fibra=[ 'cf= ',num2str(cf,5)];
mesofase=[ 'cmf= ',num2str(cmf,5)];
rho1=[ '$\rho_1$','=',num2str(cmf/cm,5),];
rho2=[ '$\rho_2$','=',num2str(cf/cm,5)];
angulo=['$\theta=$ ' num2str(theta(ii)*180/pi,3),'$^o$'];
Modulo=[ '|w2|= ',num2str(abs(w2),5)];
volmesofase=[ '$V_2= $',num2str(V2,3)];
disp(Modulo)
disp(angulo)
disp(matriz)
disp(fibra)
disp(mesofase)


Tabla(:,:)=[CENo0(:,1), CENo1(:,1),CENo2(:,1),CENo3(:,1),CENo4(:,1),CENo5(:,1),...
             CENo6(:,1), CENoG(:,1)];

 disp('Vf  \kappa_11  n=0,         1            2                3             4              5           6             Om ')    
 disp(' -------------------------------------------')
 disp(num2str([lam' Tabla(:,:)],6))
 
 disp(' -------------------------------------------')
 Tabla(:,:)=[CENo0(:,2), CENo1(:,2),CENo2(:,2),CENo3(:,2),CENo4(:,2),CENo5(:,2),...
      CENo6(:,2), CENoG(:,2)];
  
 disp('Vf \kappa_22  n=0,     1            2                3             4              5           6             Om ')  
 disp(' -------------------------------------------')
 disp(num2str([lam' Tabla(:,:)],6))
 
 Tabla(:,:)=[CENo0(:,3), CENo1(:,3),CENo2(:,3),CENo3(:,3),CENo4(:,3),CENo5(:,3),...
             CENo6(:,3), CENoG(:,3)];
 disp(' -------------------------------------------')
disp('Vf  \kappa_12  n=0,     1            2                3             4              5           6           Om ') 
 disp(' -------------------------------------------')
 disp(num2str([lam' Tabla(:,:)],6))
  
end
 
colores = [[1, 0, 0]; [0, 0, 1]; [1, 0, 1];  [0, 1, 0]; ...
           [0, 0.8, 0.8]; [.8,0.5,.8]; [0, 0, 0]; [0.5,0.5,0.5]; [0.5,0.5,1]; ...
           [0.8,0.8,0]; [1, 0.4, 0]; [.8,0,0]; [0,.8,.8]];
       
estilos_linea = {'-.', '--', '-', ':', '-.', '--', ':', '-'};  % Diferentes estilos de línea
%marcadores = {'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o'};        % Diferentes marcadores
lineawidth={2.5, 2.5, 2.5, 2.0, 1.5, 2.5, 2.5, 2.5};
letras = ['a)';'b)';'c)';'d)';'e)';'f)';'g)';'h)';'i)'];
figure(1)
  
for tt = 1:2
    for jj = 1:3
        subplot(3, 3, 3*(tt-1)+jj)
        hold on
        
        % Lista de funciones a graficar (ajusta según tus variables reales)
        funciones = {CENoG, CENo5, CENo4, CENo3, CENo2, CENo1, CENo0};
        
        for kk = 1:length(funciones)
            plot(lam(1:Lon), funciones{kk}(1:Lon,jj,tt), ...
                 'LineStyle', estilos_linea{kk}, ...
                 'Color', colores(kk,:), ...
                 'LineWidth', lineawidth{kk}, ...
                 'MarkerEdgeColor', colores(kk,:), ...
                 'MarkerFaceColor', 'w', ...
                 'MarkerSize', 4)
        end
        
        hold off

        ylabel(textoY(jj,:), 'fontsize', 18, 'interpreter', 'latex')
%         if tt == 2
%             xlabel('Fiber volume fraction ($V_3$)', 'fontsize', 12, 'interpreter', 'latex')
%         end
        xlim([0, lam(Lon)])
        xticks([0:.1:lam(Lon)]);
        text(.02, CENoG(Lon,jj,tt), [letras(3*(tt-1)+jj,:)], 'fontsize', 14, 'interpreter', 'latex')
        grid on
        grid minor
        box on
        title(['$V_2=$', num2str(me(tt),1)], 'fontsize', 12, 'interpreter', 'latex')
        
        if jj == 1 && tt == 1
            legend('$O_{30}$', '$O_5$', '$O_4$', '$O_3$', '$O_2$', '$O_1$', '$O_0$', ...
                   'fontsize', 14, 'NumColumns',7, 'interpreter', 'latex')
            grid on
            text(.01, 1.8, ['$\rho_1=$', num2str(xas(1),3)], 'fontsize', 14, 'interpreter', 'latex')
            text(.01, 1.6, ['$\rho_2=$', num2str(xas(2),3)], 'fontsize', 14, 'interpreter', 'latex')
            text(.01, 1.4, angulo, 'fontsize', 12, 'interpreter', 'latex')
        end
    end
end


textoY=['$\hat{\kappa}_{11}^{3f}/\hat{\kappa}_{11}^{2f}$';'$\hat{\kappa}_{22}^{3f}/\hat{\kappa}_{22}^{2f}$';...
       '$\hat{\kappa}_{12}^{3f}/\hat{\kappa}_{12}^{2f}$'];


    for jj = 1:3
        subplot(3, 3, jj+6)
        hold on
        
        % Lista de funciones a graficar (ajusta según tus variables reales)
        funciones = {CENoG, CENo5, CENo4, CENo3, CENo2, CENo1, CENo0};
        
        for kk = 1:length(funciones)
            plot(lam(1:Lon), funciones{kk}(1:Lon,jj,2)./funciones{kk}(1:Lon,jj,1), ...
                 'LineStyle', estilos_linea{kk}, ...
                 'Color', colores(kk,:), ...
                 'LineWidth', lineawidth{kk}, ...
                 'MarkerEdgeColor', colores(kk,:), ...
                 'MarkerFaceColor', 'w', ...
                 'MarkerSize', 4)
        end
        
        hold off

        ylabel(textoY(jj,:), 'fontsize', 18, 'interpreter', 'latex')
        if tt == 2
            xlabel('Fiber volume fraction ($V_3$)', 'fontsize', 12, 'interpreter', 'latex')
        end
        xlim([0, lam(Lon)])
        xticks([0:.1:lam(Lon)]);
        text(.02, CENoG(Lon,jj,2)/CENoG(Lon,jj,1), [letras(jj+6,:)], 'fontsize', 14, 'interpreter', 'latex')
        grid on
        grid minor
        box on
       
        
%         if jj == 1 
%             legend('$O_{30}$', '$O_5$', '$O_4$', '$O_3$', '$O_2$', '$O_1$', '$O_0$', ...
%                    'fontsize', 14, 'NumColumns',7, 'interpreter', 'latex')
%             grid on
%             text(.12,1.3,['$\rho_1=$',num2str(xas(1),3)],'fontsize',14,'interpreter','latex')
%             text(.12,1.2,['$\rho_2=$',num2str(xas(2),3)],'fontsize',14,'interpreter','latex')
%             text(.12,1.1,angulo, 'fontsize',12,'interpreter','latex')
%             text(.12,1,volmesofase,'fontsize',12, 'fontsize',12,'interpreter','latex')
%         end
    end



%%% orden del sistema No=1
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
W31=12*sqrt(21)*X7*r1^10*L10;  W32=66*sqrt(35)*X7*r1^12*L12;   W33=1716*X7*r1^14*L14;           W34=715*sqrt(63)*X7*r1^16*L16;     
W41=55*sqrt(3)*X9*r1^12*L12;  W42=429*sqrt(5)*X9*r1^14*L14;   W43=715*sqrt(63)*X9*r1^16*L16;   W44=24310*X9*r1^18*L18;  

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
X1=((1-x1)*(x1+x2)*(V2+V3)+(1+x1)*(x1-x2)*(V3))/...
             ((1+x1)*(x1+x2)*(V2+V3)+(1-x1)*(x1-x2)*(V3));
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
X1=((1-x1)*(x1+x2)*(V2+V3)+(1+x1)*(x1-x2)*(V3))/...
             ((1+x1)*(x1+x2)*(V2+V3)+(1-x1)*(x1-x2)*(V3));
S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);S16=S1(16);S18=S1(18);
S20=S1(20);S22=S1(22);S24=S1(24);S26=S1(26);

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
                    
      % Input parameters:
    % m - order of approximation (m=1,2,3,...)
    % omega1, omega2 - complex periods (omega1=1, omega2=r*exp(i*theta))
    % R2 - radius of central fiber
    % t - coating thickness
    % Sk - complex coefficients (array)
    % delta1, delta2 - parameters
    % kappa1, kappa2, kappa3 - conductivity coefficients
    

    V2plusV3=V2+V3;
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

