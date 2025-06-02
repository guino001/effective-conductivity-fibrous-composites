clearvars
%load S10_45_75_90_delta1_2; 
load Sk_Tk_30_mas_5_675_hasta_90_ENoP
ro1=[990.5, 99.498, 10.41953]; %%x1=cmf/cm;
%% 60
delta_1=delta1;
delta_2=delta2;
V2masV3=[0.1, 0.2, 0.3, 0.4, 0.5, 0.55]; 
xx=[ 1.0, 0.01, 0.0001];

lam=V2masV3;
ang=[3:1:14];
rm=min(sin(theta/2),0.5);
percolacion =pi*rm.^2./sin(theta);
angulos=(theta*180/pi);
%%% propiedades materiales
cf=1/100;
cm=1;
x2=cf/cm;
for ii=1:length(ang)
 S1= Sk(:,ang(ii));
 Periodos2=(cos(theta(ang(ii)))+1i*sin(theta(ang(ii))));
 w1=1;
 w2=Periodos2;
delta1=delta_1(ang(ii));
delta2=delta_2(ang(ii)); 

 %***********************************************************
% Metodo de Homogenizacion Asintotica
%*************************************************************
V = abs(w1)*abs(w2)*sin(theta(ang(ii))); % *** Volumen de la celda periodica *******
H1=(conj(delta1)*conj(w2)-conj(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
H2=((delta1)*conj(w2)-(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
% *************************************************************************

for jj=1:length(V2masV3)
   for kk=1:length(xx) %  t/R1
    lam3(jj,kk)=V2masV3(jj)/(1+xx(kk)).^2; % volumen de fibra interna
    lam2(jj,kk)=V2masV3(jj)-lam3(jj,kk); % volumen del anillo
    V2=lam2(jj,kk);
    V3=lam3(jj,kk);
    r1=sqrt(V*V2masV3(jj)/pi);  % radio fibra mas mesofase
    r2=sqrt(V*V3/pi);  % radio fibra mas fibra interior
    t=r1-r2;
   % *** Radio de la fibra más la mesophase ************************************
      for rr=1:length(ro1) 
      cmf=ro1(rr);
      x1=cmf/cm;
      X1=((1-x1)*(x1+x2)*(V2+V3)+(1+x1)*(x1-x2)*(V3))/...
             ((1+x1)*(x1+x2)*(V2+V3)+(1-x1)*(x1-x2)*(V3));
      J=eye(2)+X1*r1^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)];
 
       CENo(jj,:,kk,rr,ii) =ahm3f_Paralelogramo_Cortas_2025_No5(x1,x2,V3,V2,S1,X1,r1,J);
       %CENo1(jj,:,kk,rr) = ahm3f_elastico_p13_p23_Paralelogramo_Agosto(cm,cf,cmf,V3,V2,V,S1,H1,H2,31);
       %CENo(jj,:,kk,rr,ii)=conductividad_efectiva_orden_m(cm,cf,cmf,V3,V2,S1,r1,J,30);
%        CELu(jj,kk,rr)=LU1996(cm,cf,cmf,V3,V2);
%        CENo0(jj,:,kk,rr) =coeficientes_efectivos(x1,x2,V3,V2,X1,J);
end
end
end
end
%%% para comparar formula corta de Lu con O_0 nuestra


 textoxx=[];
for rr=1:length(ro1)
for kk=1:length(xx)
Tabla1=[V2masV3'];Tabla2=[V2masV3']; Tabla3=[V2masV3'];
for ii=1:length(ang)
       Tabla1=[Tabla1,CENo(:,1,kk,rr,ii)];
       Tabla2=[Tabla2,CENo(:,2,kk,rr,ii)];
       Tabla3=[Tabla3,CENo(:,3,kk,rr,ii)];
textoxx=['$t/R_1=$', num2str(xx(kk)), ',', '$\rho_1=$', num2str(ro1(rr)) ];
end
disp('--------------------------------------------------')
disp(textoxx)
disp('--------------------------------------------------')         
disp(['Comparacion entre AHM ', '$O_6$', 'y Lu 1996', '$\kappa_{11}/\kappa_1$'])
disp ('Vf       40          45        50        55       60       65     67.5       70         75         80         85         90')
 disp(num2str([ Tabla1],3))
 disp(['Comparacion entre AHM', '$O_6$', 'y Lu 1996', '$\kappa_{22}/\kappa_1$'])
 disp ('Vf       40          45        50        55       60       65     67.5       70         75         80         85        90')
 disp(num2str([ Tabla2],3))
 disp(['Comparacion entre AHM', '$O_6$', 'y Lu 1996', '$\kappa_{12}/\kappa_1$'])
 disp ('Vf       40          45        50        55       60       65     67.5       70         75          80         85        90')
 disp(num2str([ Tabla3],3))
end
end

%  
angulos=theta(ang)*180/pi;
matriz=[ 'cm= ',num2str(cm,5)];
fibra=[ 'cf= ',num2str(cf,5)];
rho2=[ '$\rho_2$','=',num2str(cf/cm,5)];
colors = lines(length(V2masV3));  % 'lines' genera colores distintos, también puedes usar 'jet', 'hsv', etc.
ejey = ['$\hat{\kappa}_{11}/\kappa_1$'; '$\hat{\kappa}_{22}/\kappa_1$'; '$\hat{\kappa}_{12}/\kappa_1$'];
colores=[[0,0,1];[.6,.6,0];[0,1,0];[1,0,0];[0,1,1];[1,0,1];[0,0,0];[0.5,0.5,0.5];[0.5,0.5,1];[0.8,0.8,0];[0.1,0.8,0];...
         [.8,0,0];[0,.8,.8];[.1,.1,0];[0,.1,0];[.1,0,0]; [0.9, 0.7, 0.5 ]; [0.75, 0.5, 0.5];[0.25, 0.25, 0.5];...
          [0.75, 0.75, 0];[0.5, 0, 0.5]];
 incisos=['a)'; 'b)'; 'c)'; 'd)'; 'e)';'f)'];
 Lu96Ram = Lu1996Tabla3();
 xram=randi([40, 90])*[1,1,1,1,1]'; %%% xram=46
 
 for rr=1:1:length(1)
 figure(rr)
for kk=1:length(xx)

%set(gcf, 'Position', [100 100 1200 400]); % Ajusta el tamaño de la figura   
   for cc = 1:3
    subplot(3,3, 3*(kk-1)+cc)
    cm=1; cmf=ro1(rr);
    rho1=[ '$\rho_1$','=',num2str(cmf/cm,5)];
    t_overR1=[['$t/R_2$='], num2str(xx(kk),3)];
    hold on
    
    % Prealocación manual para handles y leyenda
    h = zeros(length(V2masV3), 1);
    legendEntries = cell(length(V2masV3), 1);
    
    for jj = 1:length(V2masV3)
    % Seleccionar estilo basado en el índice kk
    switch mod(jj,6)
        case 1
            linestyle = '-';  linewidth = 1.5;
        case 2
            linestyle = '--'; linewidth = 1.5;
        case 3
            linestyle = ':';  linewidth = 2.0;
        case 4
            linestyle = '-.'; linewidth = 1.5;
        case 5
            linestyle = '-';  linewidth = 2.5;
        case 6
            linestyle = '--'; linewidth = 2.0;
%         case 7
%             linestyle = ':';  linewidth = 2.5;
%         case 8
%             linestyle = '-.';  linewidth = 2.5;
    end
    
    CENo_z = reshape(CENo(jj, cc, kk, rr, :), [], 1);   
    h(jj) = plot(angulos,(CENo_z), ...
                'LineStyle', linestyle, ...
                'LineWidth', linewidth, ...
                'Color', colores(jj,:));
    legendEntries{jj} = ['$V_f = $' num2str(V2masV3(jj),3)];
    end
    
      if  (cc==1 | cc==2)
      
   h(jj+1)= plot(46*[1,1,1,1,1]',Lu96Ram(:,3*(kk-1)+1,rr) ,'o','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','g',...
                       'MarkerSize',4);
    
    h(jj+2)= plot(90*[1,1,1,1,1]',Lu96Ram(:,3*(kk-1)+2,rr) ,'s','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','g',...
                       'MarkerSize',4);
              
    h(jj+3)= plot(60*[1,1,1,1,1]',Lu96Ram(:,3*(kk-1)+3,rr) ,'d','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','g',...
                       'MarkerSize',4);
      end
      legendEntries{7} = ['random']; legendEntries{8} = ['square']; legendEntries{9} = ['hexagonal'];
    hold off
    title( [rho1,', ',rho2,', ',t_overR1 ],'fontsize', 12, 'interpreter', 'latex')
    % Cuadrícula y marco
    grid on
    box on
    set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.7)
   
    
    % Etiquetas de ejes
     ylabel(ejey(cc,:), 'fontsize', 14, 'interpreter', 'latex', 'FontWeight', 'bold');
    xlabel(['$\theta^o$'], 'fontsize', 10, 'interpreter', 'latex', 'FontWeight', 'bold');
    
    if cc==1 & kk==1   % Leyenda
    legend(h, legendEntries,'FontSize', 12, 'Location', 'best',...
               'interpreter', 'latex', 'FontWeight', 'bold','NumColumns',9);
    end
    
    end
    
    % Anotaciones de texto
    %text_pos_ang = min(ang) + 0.1*(max(xx)-min(xx));
    
    % Textos fijos
    %title( [rho1,', ',rho2,', ',  angulo],'fontsize', 10, 'interpreter', 'latex')
    %title( [angulo],'fontsize', 12, 'interpreter', 'latex', 'FontWeight', 'bold');
%     text(text_pos_x, CENo(1,jj,2,ii), rho1,'fontsize', 14, 'interpreter', 'latex')
%     text(text_pos_x, CENo(1,jj,5,ii), rho2, 'fontsize', 14, 'interpreter', 'latex')
%     text(text_pos_x, CENo(1,jj,7,ii), angulo, 'fontsize', 14, 'interpreter', 'latex')
%     text(.08, max(ylim), incisos(jj,:), 'fontsize', 12, 'interpreter', 'latex')
    % Ajuste de ejes
    xlim([min(angulos) max(angulos)])
    ylim ([min(min(CENo(:, cc, kk, rr, :)))-0.01,  max(max(CENo(:, cc, kk, rr, :)))+0.01])
    
    %%%%%%%%%%%%%
    end
 end

   
%%% FUNCIONES AUXILIARES
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

RR1 =coeficientes_efectivos(x1,x2,V3,V2,X1,Z);

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

RR2 =coeficientes_efectivos(x1,x2,V3,V2,X1,Z);

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

RR3 =coeficientes_efectivos(x1,x2,V3,V2,X1,Z);

end


%%% orden del sistema No=4
function RR4 =ahm3f_Paralelogramo_Cortas_2025_No4(x1,x2,V3,V2,S1,X1,r1,J)

S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);S16=S1(16);S18=S1(18);
S20=S1(20);S22=S1(22);

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

RR4 =coeficientes_efectivos(x1,x2,V3,V2,X1,Z);

end


%%% orden del sistema No=5
function RR5 =ahm3f_Paralelogramo_Cortas_2025_No5(x1,x2,V3,V2,S1,X1,r1,J)

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

RR5 =coeficientes_efectivos(x1,x2,V3,V2,X1,Z);

end

%%% orden del sistema No=6
function RR5 =ahm3f_Paralelogramo_Cortas_2025_No6(x1,x2,V3,V2,S1,X1,r1,J)

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

RR5 =coeficientes_efectivos(x1,x2,V3,V2,X1,Z);

end



function [kappa_eff] = conductividad_efectiva_orden_m(kappa1, kappa3, kappa2,V3,V2,Sk,R1, J1,m)
                    
      % Input parameters:
    % m - order of approximation (m=1,2,3,...)
    % omega1, omega2 - complex periods (omega1=1, omega2=r*exp(i*theta))
    % R2 - radius of central fiber
    % t - coating thickness
    % Sk - complex coefficients (array)
    % delta1, delta2 - parameters
    % kappa1, kappa2, kappa3 - conductivity coefficients
    
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
    kappa_eff = coeficientes_efectivos(rho1,rho2,V3,V2,chi(1),Z);
end


%**************************************************************************
% FORMULA DE LOS COEFICIENTES EFECTIVOS (PROBLEMAS Antiplano)
%**************************************************************************
function RR =coeficientes_efectivos(x1,x2,V3,V2,X1,Z)
F=(1-x1+2*V3*(x1-x2)*x1/(V2*(x1+x2)+2*V3*x1));
G=1/(V2*(x1+x2)+2*V3*x1);
DtZ=det(Z);
c55=(1-V2-V3+x1*V2+x2*V3-(V2+V3)*F*((X1+1)*Z(2,2)-DtZ)/(DtZ)-G*(x1-x2)^2*V2*V3);
                    c45=((V2+V3)*F*((X1+1)*Z(2,1))/(DtZ));
                    c54=((V2+V3)*F*((X1+1)*Z(1,2))/(DtZ));
c44=(1-V2-V3+x1*V2+x2*V3-(V2+V3)*F*((X1+1)*Z(1,1)-DtZ)/(DtZ)-G*(x1-x2)^2*V2*V3);
      
RR=[c55,c44,c54,c45];
end

%%%% tabla 3 pag 4402
 function Lu96Ram=Lu1996Tabla3
                 %%%     t/R_=1     t/R_=0.01        t/R_=0.0001
                %%% Ra   Sa   Ha   Ra   Sa   Ha   Ra   Sa   Ha
 Lu96Ram(:,:,1)= [ 1.23 1.22 1.22 1.18 1.18 1.18 0.848 0.851 0.851 
                   1.54 1.50 1.50 1.41 1.39 1.39 0.709 0.723 0.723 
                   1.98 1.86 1.85 1.71 1.65 1.65 0.580 0.611 0.612 
                   2.62 2.34 2.33 2.11 1.98 1.97 0.455 0.512 0.514 
                   3.60 3.07 2.99 2.66 2.41 2.38 0.332 0.424 0.427];

Lu96Ram(:,:,2)=  [ 1.22 1.21 1.21 1.00 1.00 1.00 0.818 0.825 0.825 
                   1.52 1.48 1.48 1.00 1.00 1.00 0.652 0.677 0.678 
                   1.93 1.82 1.82 1.00 1.00 1.00 0.493 0.552 0.552 
                   2.52 2.28 2.26 1.00 1.00 1.00 0.336 0.442 0.445 
                   3.41 2.94 2.88 1.00 1.00 1.00 0.172 0.343 0.350];

Lu96Ram(:,:,3)= [1.16 1.16 1.16 0.849 0.853 0.853 0.815 0.822 0.822 
                 1.35 1.34 1.34 0.712 0.725 0.725 0.645 0.673 0.673 
                 1.60 1.56 1.56 0.584 0.614 0.615 0.483 0.545 0.546 
                 1.91 1.82 1.82 0.460 0.516 0.517 0.320 0.435 0.437 
                 2.31 2.16 2.14 0.338 0.427 0.430 0.150 0.335 0.343 ];
 end



%**************************************************************************
function LU =LU1996(cm,cf,cmf,V3,V2) %% es 1995
Bt13=(cm-cmf)/(cm+cmf);
Bt23=(cf-cmf)/(cf+cmf);

LU=cm*(1+2*(V2+V3)*(Bt23-Bt13*(V2+V3)/V3)/...
                   ((V2+V3)/V3 -Bt13*Bt23 - (V2+V3)*(Bt23-Bt13*(V2+V3)/V3)));
end

