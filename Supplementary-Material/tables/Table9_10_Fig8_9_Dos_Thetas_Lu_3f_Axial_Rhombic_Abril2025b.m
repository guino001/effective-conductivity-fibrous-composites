clearvars
load S10_45_75_90_delta1_2; 
ro1=[990.5, 99.498, 10.41953]; %%x1=cmf/cm;
%% 60
ii=7;
V2masV3=[0.1, 0.2, 0.3, .4, .5,.9, .905, .906899]; 
xx=[ 1.0, 0.1, 0.01, 0.001];

lam=V2masV3;
 S1=S10_45_75_90_delta1_2(:,ii);
 delta1=2*ZZ(ii);
 Periodos2=(cos(theta(ii))+1i*sin(theta(ii)));
 w1=1;
 w2=Periodos2;
rm=min(sin(theta(ii)/2),0.5);
percolacion =pi*rm.^2./sin(theta(ii));
delta1=delta1;
delta2=delta1.*w2-2*pi*1i; 

 %***********************************************************
% Metodo de Homogenizacion Asintotica
%*************************************************************
V = abs(w1)*abs(w2)*sin(theta(ii)); % *** Volumen de la celda periodica *******
H1=(conj(delta1)*conj(w2)-conj(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
H2=((delta1)*conj(w2)-(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
% *************************************************************************
%%% propiedades materiales
cf=1/100;
cm=1;
x2=cf/cm;

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
 
       CENo(jj,:,kk,rr) =ahm3f_Paralelogramo_Cortas_2025_No6(x1,x2,V3,V2,S1,X1,r1,J);
       %CENo1(jj,:,kk,rr) = ahm3f_elastico_p13_p23_Paralelogramo_Agosto(cm,cf,cmf,V3,V2,V,S1,H1,H2,31);
       CENo1(jj,:,kk,rr)=conductividad_efectiva_orden_m(cm,cf,cmf,V3,V2,S1,r1,J,30);
       CELu(jj,kk,rr)=LU1996(cm,cf,cmf,V3,V2);
       CENo0(jj,:,kk,rr) =coeficientes_efectivos(x1,x2,V3,V2,J);
end
end
end

% %%% para comparar formula corta de Lu con O_0 nuestra
% Tabla=[];
% for kk=1:length(xx)   
% Tabla=[Tabla,CENo0(:,1,kk,1), CELu(:,kk,1)];
% %textoxx=[textoxx; num2str(xx(kk))];
% end
% disp(' Comparacion entre AHM O_0 y Lu')
%  disp(num2str([lam' Tabla],6))
% %  
%% %Tabla
 Lu95hex=Lu1995Tabla4;
disp('Comparacion Lu1995 con spring AHM')
% Configuración de formato
decimales_lam = 6;   % Decimales para lambda
decimales_val = 6;   % Decimales para valores AHM/Lu
ancho_col = 15;      % Ancho de cada columna

num_bi = length(xx); % Debe ser par (6 en este caso)

for rr = 1:length(ro1)
    fprintf('\rrho_1= %.4f\n', ro1(rr))
    
    % Primera mitad (Bi1-Bi3)
    %fprintf('\nPrimera Parte (Bi1-Bi3):\n');
    % Primero: Imprimir los valores de Bi (centrados sobre AHM/Lu)
    disp('-----------------------------------------------------------------------------------------------------------')
    fprintf(['%', num2str(ancho_col), 's'], ' '); % Espacio para Lambda
    for kk = 1:num_bi/2
        fprintf(['%', num2str(ancho_col), '.1e%', num2str(ancho_col), '.1e'], xx(kk), xx(kk));
    end
    fprintf('\n');
    % Segundo: Imprimir encabezados (Lambda, AHM1, Lu1, AHM2, Lu2, ...)
    fprintf(['%', num2str(ancho_col), 's'], 'V_3');
    for kk = 1:num_bi/2
        fprintf(['%', num2str(ancho_col), 's%', num2str(ancho_col), 's'], ...
                sprintf('AHM', kk), sprintf('Lu', kk));
    end
    fprintf('\n%s\n', repmat('-', 1, (1 + num_bi)*ancho_col));
    disp('        ')
    % Tercero: Imprimir los datos
    for ii = 1:length(lam)
        fprintf(['%', num2str(ancho_col), '.', num2str(decimales_lam), 'f'], lam(ii));
        for kk = 1:num_bi/2
            fprintf(['%', num2str(ancho_col), '.', num2str(decimales_val), 'f'], CENo(ii,1,kk,rr));
            fprintf(['%', num2str(ancho_col), '.', num2str(decimales_val), 'f'], Lu95hex(ii,kk,rr));
        end
        fprintf('\n');
    end
    
    % Segunda mitad (Bi4-Bi6)
    %fprintf('\nSegunda Parte (Bi4-Bi6):\n');
    disp('-----------------------------------------------------------------------------------------------------------')
    % Primero: Imprimir los valores de Bi (centrados sobre AHM/Lu)
    fprintf(['%', num2str(ancho_col), 's'], ' '); % Espacio para Lambda
    for kk = num_bi/2+1:num_bi
        fprintf(['%', num2str(ancho_col), '.1e%', num2str(ancho_col), '.1e'], xx(kk), xx(kk));
    end
    fprintf('\n');
    
    % Segundo: Imprimir encabezados (Lambda, AHM4, Lu4, AHM5, Lu5, ...)
    fprintf(['%', num2str(ancho_col), 's'], 'Lambda');
    for kk = num_bi/2+1:num_bi
        fprintf(['%', num2str(ancho_col), 's%', num2str(ancho_col), 's'], ...
                sprintf('AHM', kk), sprintf('Lu', kk));
    end
    fprintf('\n%s\n', repmat('-', 1, (1 + num_bi)*ancho_col));
    
    % Tercero: Imprimir los datos
    for ii = 1:length(lam)
        fprintf(['%', num2str(ancho_col), '.', num2str(decimales_lam), 'f'], lam(ii));
        for kk = num_bi/2+1:num_bi
            fprintf(['%', num2str(ancho_col), '.', num2str(decimales_val), 'f'], CENo(ii,1,kk,rr));
            fprintf(['%', num2str(ancho_col), '.', num2str(decimales_val), 'f'], Lu95hex(ii,kk,rr));
        end
        fprintf('\n');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%
disp('%%%%%%%%%%%%%%%%%%%%% ')
disp('Comparacion con la Exacta JuanK 2011 para Vf alto')
% Configuración de formato
decimales_lam = 6;   % Decimales para lambda
decimales_val = 6;   % Decimales para valores AHM/Lu
ancho_col = 15;      % Ancho de cada columna
% Bi = [1.e12, 1, 0.1, 0.01, 0.001, 0.0001]; % Valores de Bi
num_bi = length(xx); % Debe ser par (6 en este caso)

for rr = 1:length(ro1)
    fprintf('\rrho= %.4f\n', ro1(rr))
     disp('-----------------------------------------------------------------------------------------------------------')
    % Primera mitad (Bi1-Bi3)
    %fprintf('\nPrimera Parte (Bi1-Bi3):\n');
    % Primero: Imprimir los valores de Bi (centrados sobre AHM/Lu)
    fprintf(['%', num2str(ancho_col), 's'], ' '); % Espacio para Lambda
    for kk = 1:num_bi/2
        fprintf(['%', num2str(ancho_col), '.1e%', num2str(ancho_col), '.1e'], xx(kk), xx(kk));
    end
    fprintf('\n');
    % Segundo: Imprimir encabezados (Lambda, AHM1, Lu1, AHM2, Lu2, ...)
    fprintf(['%', num2str(ancho_col), 's'], 'V_3');
    for kk = 1:num_bi/2
        fprintf(['%', num2str(ancho_col), 's%', num2str(ancho_col), 's'], ...
                sprintf('AHM', kk), sprintf('Lu', kk));
    end
    fprintf('\n%s\n', repmat('-', 1, (1 + num_bi)*ancho_col));
    disp('        ')
    % Tercero: Imprimir los datos
    for ii = 6:length(lam)
        fprintf(['%', num2str(ancho_col), '.', num2str(decimales_lam), 'f'], lam(ii));
        for kk = 1:num_bi/2
            fprintf(['%', num2str(ancho_col), '.', num2str(decimales_val), 'f'], CENo1(ii,1,kk,rr));
            fprintf(['%', num2str(ancho_col), '.', num2str(decimales_val), 'f'], Lu95hex(ii,kk,rr));
        end
        fprintf('\n');
    end
    
    % Segunda mitad (Bi4-Bi6)
    %fprintf('\nSegunda Parte (Bi4-Bi6):\n');
     disp('-----------------------------------------------------------------------------------------------------------')
    % Primero: Imprimir los valores de Bi (centrados sobre AHM/Lu)
    fprintf(['%', num2str(ancho_col), 's'], ' '); % Espacio para Lambda
    for kk = num_bi/2+1:num_bi
        fprintf(['%', num2str(ancho_col), '.1e%', num2str(ancho_col), '.1e'], xx(kk), xx(kk));
    end
    fprintf('\n');
    
    % Segundo: Imprimir encabezados (Lambda, AHM4, Lu4, AHM5, Lu5, ...)
    fprintf(['%', num2str(ancho_col), 's'], 'V_3');
    for kk = num_bi/2+1:num_bi
        fprintf(['%', num2str(ancho_col), 's%', num2str(ancho_col), 's'], ...
                sprintf('AHM', kk), sprintf('Lu', kk));
    end
    fprintf('\n%s\n', repmat('-', 1, (1 + num_bi)*ancho_col));
    
    % Tercero: Imprimir los datos
    for ii = 6:length(lam)
        fprintf(['%', num2str(ancho_col), '.', num2str(decimales_lam), 'f'], lam(ii));
        for kk = num_bi/2+1:num_bi
            fprintf(['%', num2str(ancho_col), '.', num2str(decimales_val), 'f'], CENo1(ii,1,kk,rr));
            fprintf(['%', num2str(ancho_col), '.', num2str(decimales_val), 'f'], Lu95hex(ii,kk,rr));
        end
        fprintf('\n');
    end
end
%%%%% Comparacion con dos celdas diferentes y los tres constrastes de Lu
%%%%%%%%%%%%%%%%%%%%%%%%% primer angulo
load S10_45_75_90_delta1_2; 
%ro1=[990.5, 99.498, 10.41953]; %%x1=cmf/cm;
%%%%% primera celda 
ii=5; % 45
xx=[0:.0001: 0.01, .02:.02:.2, .21:.1:.9, 1]; 
V2masV3=[.65,.64,.62, .6:-.1:.2];
lam=V2masV3;
lam1(1,:)=lam;
cm=1;
%cf=100;
x2=cf/cm;
S1=S10_45_75_90_delta1_2(:,ii);
delta1=2*ZZ(ii);
Periodos2=(cos(theta(ii))+1i*sin(theta(ii)));
w1=1;
w2=Periodos2;
rm=min(sin(theta(ii)/2),0.5);
percolacion =pi*rm.^2./sin(theta(ii));
delta1=delta1;
delta2=delta1.*w2-2*pi*1i; 
 %***********************************************************
% Metodo de Homogenizacion Asintotica
%*************************************************************
V = abs(w1)*abs(w2)*sin(theta(ii)); % *** Volumen de la celda periodica *******
angulo(1,:)=['$\theta$ =' num2str(theta(ii)*180/pi,3),'$^o$'];
H1=(conj(delta1)*conj(w2)-conj(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
H2=((delta1)*conj(w2)-(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);

for jj=1:length(V2masV3)
   for kk=1:length(xx)
    lam3(jj,kk)=V2masV3(jj)/((1+xx(kk)).^2); % volumen de fibra interna
    lam2(jj,kk)=V2masV3(jj)-lam3(jj,kk); % volumen del anillo
    V2=lam2(jj,kk);
    V3=lam3(jj,kk);
    r1=sqrt(V*V2masV3(jj)/pi);  % radio fibra mas mesofase
   % *** Radio de la fibra más la mesophase ************************************
       for rr=1:length(ro1)
        cmf=ro1(rr);
        x1=cmf/cm;
        X1=((1-x1)*(x1+x2)*(V2+V3)+(1+x1)*(x1-x2)*(V3))/...
             ((1+x1)*(x1+x2)*(V2+V3)+(1-x1)*(x1-x2)*(V3));
        J=eye(2)+X1*r1^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)];
       
       CENo(jj,:,kk,rr,1) =ahm3f_Paralelogramo_Cortas_2025_No6(x1,x2,V3,V2,S1,X1,r1,J);
       %CENo1(jj,:,kk,rr,1) = ahm3f_elastico_p13_p23_Paralelogramo_Agosto(cm,cf,cmf,V3,V2,V,S1,H1,H2,31);
       %CENo(jj,:,kk,rr,1)=conductividad_efectiva_orden_m(cm, cf, cmf,V3,V2,V,S1,H1,H2,15)
        end
   end
end
 %%%%% Segunda  celda 
ii=9; % 75
lam=[.81,.8:-.1:.2];
%lam=[.75, .72 .7:-.1:.2];
lam1(2,:)=lam;
S1=S10_45_75_90_delta1_2(:,ii);
delta1=2*ZZ(ii);
Periodos2=(cos(theta(ii))+1i*sin(theta(ii)));
w1=1;
w2=Periodos2;
rm=min(sin(theta(ii)/2),0.5);
percolacion =pi*rm.^2./sin(theta(ii));
delta1=delta1;
delta2=delta1.*w2-2*pi*1i; 
 %***********************************************************
% Metodo de Homogenizacion Asintotica
%*************************************************************
V = abs(w1)*abs(w2)*sin(theta(ii)); % *** Volumen de la celda periodica *******
angulo(2,:)=['$\theta$ =' num2str(theta(ii)*180/pi,3),'$^o$'];
H1=(conj(delta1)*conj(w2)-conj(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
H2=((delta1)*conj(w2)-(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
% *************************************************************************

for jj=1:length(V2masV3)
   for kk=1:length(xx)
    lam3(jj,kk)=V2masV3(jj)/((1+xx(kk)).^2); % volumen de fibra interna
    lam2(jj,kk)=V2masV3(jj)-lam3(jj,kk); % volumen del anillo
    V2=lam2(jj,kk);
    V3=lam3(jj,kk);
    r1=sqrt(V*V2masV3(jj)/pi);  % radio fibra mas mesofase
   % *** Radio de la fibra más la mesophase ************************************
       for rr=1:length(ro1)
        cmf=ro1(rr);
        x1=cmf/cm;
        X1=((1-x1)*(x1+x2)*(V2+V3)+(1+x1)*(x1-x2)*(V3))/...
             ((1+x1)*(x1+x2)*(V2+V3)+(1-x1)*(x1-x2)*(V3));
        J=eye(2)+X1*r1^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)];
       
       CENo(jj,:,kk,rr,2) =ahm3f_Paralelogramo_Cortas_2025_No6(x1,x2,V3,V2,S1,X1,r1,J);
       %CENo1(jj,:,kk,rr,2) = ahm3f_elastico_p13_p23_Paralelogramo_Agosto(cm,cf,cmf,V3,V2,V,S1,H1,H2,31);
       %CENo(jj,:,kk,rr,2)=conductividad_efectiva_orden_m(cm, cf, cmf,V3,V2,V,S1,H1,H2,15);
       end
end
end

matriz=[ 'cm= ',num2str(cm,5)];
fibra=[ 'cf= ',num2str(cf,5)];
% mesofase=[ 'cmf= ',num2str(cmf,5)];
% rho1=[ '$\rho_1$','=',num2str(cmf/cm,5),];
rho2=[ '$\rho_2$','=',num2str(cf/cm,5)];
% angulo=['$\theta$ =' num2str(theta(ii)*180/pi,3),'$^o$'];
Modulo=[ '|w2|= ',num2str(abs(w2),5)];
volmesofase=[ '$V_2= $',num2str(V2,3)];

% Definir una paleta de colores (ej. usando 'lines', 'parula', 'jet', o personalizada)
colors = lines(length(xx));  % 'lines' genera colores distintos, también puedes usar 'jet', 'hsv', etc.
ejey = ['$\hat{\kappa}_{11}/\kappa_1$'; '$\hat{\kappa}_{22}/\kappa_1$'; '$\hat{\kappa}_{12}/\kappa_1$'];


colores=[[0,0,1];[.6,.6,0];[0,1,0];[1,0,0];[0,1,1];[1,0,1];[0,0,0];[0.5,0.5,0.5];[0.5,0.5,1];[0.8,0.8,0];[0.1,0.8,0];...
         [.8,0,0];[0,.8,.8];[.1,.1,0];[0,.1,0];[.1,0,0]; [0.9, 0.7, 0.5 ]; [0.75, 0.5, 0.5];[0.25, 0.25, 0.5];...
          [0.75, 0.75, 0];[0.5, 0, 0.5]];
 incisos=['a)'; 'b)'; 'c)'; 'd)'; 'e)';'f)'];
for rr=1:length(1)
figure(rr)
for ii=1:2
     lam=lam1(ii,:);
%set(gcf, 'Position', [100 100 1200 400]); % Ajusta el tamaño de la figura   
   for jj = 1:3
    subplot(2,6, 3*(ii-1)+jj)
    cm=1; cmf=ro1(rr);
    rho1=[ '$\rho_1$','=',num2str(cmf/cm,5)];
    angulo=['$\theta$ =' num2str(theta(4*ii+1)*180/pi,3),'$^o$'];
    hold on
    
    % Prealocación manual para handles y leyenda
    h = zeros(length(lam), 1);
    legendEntries = cell(length(lam), 1);
    
    for kk = 1:length(lam)
    % Seleccionar estilo basado en el índice kk
    switch mod(kk,8)
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
        case 7
            linestyle = ':';  linewidth = 2.5;
        case 8
            linestyle = '-.';  linewidth = 2.5;
    end
    
    CENo_z = reshape(CENo(kk, jj, :, rr, ii), [], 1);   
    h(kk) = plot(xx, CENo_z, ...
                'LineStyle', linestyle, ...
                'LineWidth', linewidth, ...
                'Color', colores(kk,:));
    legendEntries{kk} = ['$V_f = $' num2str(lam(kk),3)];
    end
 
    hold off
    
    % Etiquetas de ejes
     ylabel(ejey(jj,:), 'fontsize', 14, 'interpreter', 'latex', 'FontWeight', 'bold');
    xlabel(['$t/R_2$'], 'fontsize', 10, 'interpreter', 'latex', 'FontWeight', 'bold');
    
    if jj==1
    % Leyenda
    legend(h, legendEntries, 'FontSize', 10, 'Location', 'best', 'interpreter', 'latex', 'FontWeight', 'bold','NumColumns',4);
    title( [rho1,', ',rho2,', ',  angulo],'fontsize', 12, 'interpreter', 'latex')
    end
    
    % Cuadrícula y marco
    grid on
    box on
    set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.7)
    
    % Anotaciones de texto
    text_pos_x = min(xx) + 0.1*(max(xx)-min(xx));
    
    % Textos fijos
    %title( [rho1,', ',rho2,', ',  angulo],'fontsize', 10, 'interpreter', 'latex')
    %title( [angulo],'fontsize', 12, 'interpreter', 'latex', 'FontWeight', 'bold');
%     text(text_pos_x, CENo(1,jj,2,ii), rho1,'fontsize', 14, 'interpreter', 'latex')
%     text(text_pos_x, CENo(1,jj,5,ii), rho2, 'fontsize', 14, 'interpreter', 'latex')
%     text(text_pos_x, CENo(1,jj,7,ii), angulo, 'fontsize', 14, 'interpreter', 'latex')
%     text(.08, max(ylim), incisos(jj,:), 'fontsize', 12, 'interpreter', 'latex')
    % Ajuste de ejes
    xlim([min(xx) max(xx)])
    %ylim auto
    
    %%%%%%%%%%%%%
    end
   end
end

for rr=1:length(1)
figure(rr)
for ii=1:2
     lam=lam1(ii,:);
%set(gcf, 'Position', [100 100 1200 400]); % Ajusta el tamaño de la figura   
   for jj = 1:3
    subplot(2,6, 3*(ii-1) + jj + 6)
    cm=1; cmf=ro1(rr);
    rho1=[ '$\rho_1$','=',num2str(cmf/cm,5)];
    angulo=['$\theta$ =' num2str(theta(4*ii+1)*180/pi,3),'$^o$'];
    hold on
    
    % Prealocación manual para handles y leyenda
    h = zeros(length(lam), 1);
    legendEntries = cell(length(lam), 1);
    
     for kk = 1:length(lam)
    % Seleccionar estilo basado en el índice kk
    switch mod(kk,8)
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
        case 7
            linestyle = ':';  linewidth = 2.5;
        case 8
            linestyle = '-.';  linewidth = 2.5;
    end
    
    CENo_z = reshape(CENo(kk, jj, :, rr, ii), [], 1);   
    h(kk) = plot(xx, CENo_z, ...
                'LineStyle', linestyle, ...
                'LineWidth', linewidth, ...
                'Color', colores(kk,:));
    legendEntries{kk} = ['$V_f = $' num2str(lam(kk),3)];
    end
    hold off
    
    % Etiquetas de ejes
     ylabel(ejey(jj,:), 'fontsize', 14, 'interpreter', 'latex', 'FontWeight', 'bold');
   xlabel(['$t/R_2$'], 'fontsize', 10, 'interpreter', 'latex', 'FontWeight', 'bold');
    
    
    if jj==1
    % Leyenda
    legend(h, legendEntries, 'FontSize', 10, 'Location', 'best', 'interpreter', 'latex', 'FontWeight', 'bold','NumColumns',4);
    title( [rho1,', ',rho2,', ',  angulo],'fontsize', 12, 'interpreter', 'latex')
    end
    
    % Cuadrícula y marco
    grid on
    box on
    set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.7)
    
    % Anotaciones de texto
    text_pos_x = min(xx) + 0.1*(max(xx)-min(xx));
    
    % Textos fijos
    %title( [rho1,', ',rho2,', ',  angulo],'fontsize', 10, 'interpreter', 'latex')
    %title( [ angulo],'fontsize', 12, 'interpreter', 'latex', 'FontWeight', 'bold');
%     text(text_pos_x, CENo(1,jj,2,ii), rho1,'fontsize', 14, 'interpreter', 'latex')
%     text(text_pos_x, CENo(1,jj,5,ii), rho2, 'fontsize', 14, 'interpreter', 'latex')
%     text(text_pos_x, CENo(1,jj,7,ii), angulo, 'fontsize', 14, 'interpreter', 'latex')
%     text(.08, max(ylim), incisos(jj,:), 'fontsize', 12, 'interpreter', 'latex')
    % Ajuste de ejes
    xlim([0 .002])
    %ylim auto
    end
   end
end

%%%%%%%%%%%%%
%%%% Caso donde rho1 y rho2 son mayores que 1 y otro donde rho 1 y rho 2
%%%% son menosres que 1.
%%%%%%%%%%%%%%%%%%%%%%%%% primer angulo
clear
cm=1;
cf1=[2, 1/2];
ro1=[10, 1/10]; %%x1=cmf/cm;
x21=cf1/cm;
load S10_45_75_90_delta1_2; 
%ro1=[990.5, 99.498, 10.41953]; %%x1=cmf/cm;
%%%%% primera celda 

ii=5; % 45
xx=[0:.001: 0.01, .02:.02:.2, .21:.1:.9, 1]; 
V2masV3=[.65,.64,.62, .6:-.1:.2];
lam=V2masV3;
lam1(1,:)=lam;
cm=1;

S1=S10_45_75_90_delta1_2(:,ii);
delta1=2*ZZ(ii);
Periodos2=(cos(theta(ii))+1i*sin(theta(ii)));
w1=1;
w2=Periodos2;
rm=min(sin(theta(ii)/2),0.5);
percolacion =pi*rm.^2./sin(theta(ii));
delta1=delta1;
delta2=delta1.*w2-2*pi*1i; 
 %***********************************************************
% Metodo de Homogenizacion Asintotica
%*************************************************************
V = abs(w1)*abs(w2)*sin(theta(ii)); % *** Volumen de la celda periodica *******
angulo(1,:)=['$\theta$ =' num2str(theta(ii)*180/pi,3),'$^o$'];
H1=(conj(delta1)*conj(w2)-conj(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
H2=((delta1)*conj(w2)-(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);

for jj=1:length(V2masV3)
   for kk=1:length(xx)
    lam3(jj,kk)=V2masV3(jj)/((1+xx(kk)).^2); % volumen de fibra interna
    lam2(jj,kk)=V2masV3(jj)-lam3(jj,kk); % volumen del anillo
    V2=lam2(jj,kk);
    V3=lam3(jj,kk);
    r1=sqrt(V*V2masV3(jj)/pi);  % radio fibra mas mesofase
   % *** Radio de la fibra más la mesophase ************************************
       for rr=1:2
        cmf=ro1(rr);
         x1=cmf/cm;
         x2=x21(rr);
        X1=((1-x1)*(x1+x2)*(V2+V3)+(1+x1)*(x1-x2)*(V3))/...
             ((1+x1)*(x1+x2)*(V2+V3)+(1-x1)*(x1-x2)*(V3));
        J=eye(2)+X1*r1^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)];
       
       CENo(jj,:,kk,rr,rr,1) =ahm3f_Paralelogramo_Cortas_2025_No6(x1,x2,V3,V2,S1,X1,r1,J);
%        CENo1(jj,:,kk,rr,1) = ahm3f_elastico_p13_p23_Paralelogramo_Agosto(cm,cf,cmf,V3,V2,V,S1,H1,H2,31);
        end
   end
end
 %%%%% Segunda  celda 
ii=9; % 75
lam=[.81,.8:-.1:.2];
%lam=[.75, .72 .7:-.1:.2];
lam1(2,:)=lam;
S1=S10_45_75_90_delta1_2(:,ii);
delta1=2*ZZ(ii);
Periodos2=(cos(theta(ii))+1i*sin(theta(ii)));
w1=1;
w2=Periodos2;
rm=min(sin(theta(ii)/2),0.5);
percolacion =pi*rm.^2./sin(theta(ii));
delta1=delta1;
delta2=delta1.*w2-2*pi*1i; 
 %***********************************************************
% Metodo de Homogenizacion Asintotica
%*************************************************************
V = abs(w1)*abs(w2)*sin(theta(ii)); % *** Volumen de la celda periodica *******
angulo(2,:)=['$\theta$ =' num2str(theta(ii)*180/pi,3),'$^o$'];
H1=(conj(delta1)*conj(w2)-conj(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
H2=((delta1)*conj(w2)-(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
% *************************************************************************

for jj=1:length(V2masV3)
   for kk=1:length(xx)
    lam3(jj,kk)=V2masV3(jj)/((1+xx(kk)).^2); % volumen de fibra interna
    lam2(jj,kk)=V2masV3(jj)-lam3(jj,kk); % volumen del anillo
    V2=lam2(jj,kk);
    V3=lam3(jj,kk);
    r1=sqrt(V*V2masV3(jj)/pi);  % radio fibra mas mesofase
   % *** Radio de la fibra más la mesophase ************************************
       for rr=1:2
        cmf=ro1(rr);
        x1=cmf/cm;
         x2=x21(rr);
        X1=((1-x1)*(x1+x2)*(V2+V3)+(1+x1)*(x1-x2)*(V3))/...
             ((1+x1)*(x1+x2)*(V2+V3)+(1-x1)*(x1-x2)*(V3));
        J=eye(2)+X1*r1^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)];
       
       CENo(jj,:,kk,rr,rr,2) =ahm3f_Paralelogramo_Cortas_2025_No6(x1,x2,V3,V2,S1,X1,r1,J);
       %CENo1(jj,:,kk,rr,rr,2) = ahm3f_elastico_p13_p23_Paralelogramo_Agosto(cm,cf,cmf,V3,V2,V,S1,H1,H2,31);
       end
end
end

% Definir una paleta de colores (ej. usando 'lines', 'parula', 'jet', o personalizada)
colors = lines(length(xx));  % 'lines' genera colores distintos, también puedes usar 'jet', 'hsv', etc.
ejey = ['$\hat{\kappa}_{11}/\kappa_1$'; '$\hat{\kappa}_{22}/\kappa_1$'; '$\hat{\kappa}_{12}/\kappa_1$'];


colores=[[0,0,1];[.6,.6,0];[0,1,0];[1,0,0];[0,1,1];[1,0,1];[0,0,0];[0.5,0.5,0.5];[0.5,0.5,1];[0.8,0.8,0];[0.1,0.8,0];...
         [.8,0,0];[0,.8,.8];[.1,.1,0];[0,.1,0];[.1,0,0]; [0.9, 0.7, 0.5 ]; [0.75, 0.5, 0.5];[0.25, 0.25, 0.5];...
          [0.75, 0.75, 0];[0.5, 0, 0.5]];
 incisos=['a)'; 'b)'; 'c)'; 'd)'; 'e)';'f)'];
 figure(2)
for rr=1:2
for ii=1:2
     lam=lam1(ii,:);
%set(gcf, 'Position', [100 100 1200 400]); % Ajusta el tamaño de la figura   
   for jj = 1:3
    subplot(2,6, 3*(ii-1)+jj + 6*(rr-1) )
    cm=1; cmf=ro1(rr); cf=cf1(rr);
    rho1=[ '$\rho_1$','=',num2str(cmf/cm,5)];
    rho2=[ '$\rho_2$','=',num2str(cf/cm,5)];
    angulo=['$\theta$ =' num2str(theta(4*ii+1)*180/pi,3),'$^o$'];
    hold on
    
    % Prealocación manual para handles y leyenda
    h = zeros(length(lam), 1);
    legendEntries = cell(length(lam), 1);
    
    for kk = 1:length(lam)
    % Seleccionar estilo basado en el índice kk
    switch mod(kk,8)
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
        case 7
            linestyle = ':';  linewidth = 2.5;
        case 8
            linestyle = '-.';  linewidth = 2.5;
    end
    
    CENo_z = reshape(CENo(kk, jj, :, rr,rr, ii), [], 1);   
    h(kk) = plot(xx, CENo_z, ...
                'LineStyle', linestyle, ...
                'LineWidth', linewidth, ...
                'Color', colores(kk,:));
    legendEntries{kk} = ['$V_f = $' num2str(lam(kk),3)];
    end
 
    hold off
    
    % Etiquetas de ejes
     ylabel(ejey(jj,:), 'fontsize', 14, 'interpreter', 'latex', 'FontWeight', 'bold');
    xlabel(['$t/R_2$'], 'fontsize', 10, 'interpreter', 'latex', 'FontWeight', 'bold');
    
    if jj==1
    % Leyenda
    legend(h, legendEntries, 'FontSize', 10, 'Location', 'best', 'interpreter', 'latex', 'FontWeight', 'bold','NumColumns',4);
    title( [rho1,', ',rho2,', ',  angulo],'fontsize', 12, 'interpreter', 'latex')
    end
    
    % Cuadrícula y marco
    grid on
    box on
    set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.7)
    
    % Anotaciones de texto
    text_pos_x = min(xx) + 0.1*(max(xx)-min(xx));
    xlim([0 1])
    % Textos fijos
    %title( [rho1,', ',rho2,', ',  angulo],'fontsize', 10, 'interpreter', 'latex')
    %title( [angulo],'fontsize', 12, 'interpreter', 'latex', 'FontWeight', 'bold');
%     text(text_pos_x, CENo(1,jj,2,ii), rho1,'fontsize', 14, 'interpreter', 'latex')
%     text(text_pos_x, CENo(1,jj,5,ii), rho2, 'fontsize', 14, 'interpreter', 'latex')
%     text(text_pos_x, CENo(1,jj,7,ii), angulo, 'fontsize', 14, 'interpreter', 'latex')
%     text(.08, max(ylim), incisos(jj,:), 'fontsize', 12, 'interpreter', 'latex')
% if jj==1
%     % Ajuste de ejes
%     xlim([0 .01])
%     %ylim auto
% else
%     xlim([0 .1])
% end   
    %%%%%%%%%%%%%
    end
   end
end

%%%%%%%%%%%%%


 
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

RR4 =coeficientes_efectivos(x1,x2,V3,V2,Z);

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

RR5 =coeficientes_efectivos(x1,x2,V3,V2,Z);

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

RR5 =coeficientes_efectivos(x1,x2,V3,V2,Z);

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

%%%% tabla 2 pag 2621
 function Lu95hex=Lu1995Tabla4
 
 Lu95hex(:,:,1)= [ 1.4979	1.48693 	1.39002,	1,	0.723141;
	             2.99126	2.92195	    2.38028,	1,	0.426657;
	             5.75573	5.47952 	3.70694,	1,	0.278375;
	             9.79017	8.99004	    4.94321,	1,	0.213437;
	             24.3004	19.7183	    6.81906,	1,	0.164066;
                 50.4946	33.6087	    7.73265,	1,	0.151879;
	             90.38	    46.959, 	8.1039, 	1,	0.148832;
		            0	       0         8.34, 	    1,	0.147674];

Lu95hex(:,:,2)=  [1.47958	1.386,	    1,	 0.723336	0.677603;
                  2.87648	2.3602,	    1,	0.426992	0.35049;
	              5.30437	3.65194,	1,	0.278758	0.19188;
	               8.50552	4.84106,	1,	0.213841	0.121481;
	              17.3242	6.61433,	1,	0.164493	0.063647;
	              26.5502	7.45757,	1,	0.152313	0.046562;
	              33.107	7.78955,	1,	0.149268	0.041658;
	                38.9	7.98171,	1,	0.148111	0.03967];

Lu95hex(:,:,3)= [1.33888,	1,	0.725369,	0.677855,	0.672755;
	             2.13714,	1,	0.430491,	0.3509,	    0.342617;
                 3.08019,	1,	0.282762,	0.192343,	0.182989;
                 3.84051,	1,	0.21807,	0.121979,	0.111891;
	             4.7805,	1,	0.168952,	0.064235,	0.052056;
	             5.12111,	1,	0.156849,	0.047241,	0.032339;
	             5.22147,	1,	0.153825,	0.042386,	0.025438,
	             5.26222	1,	0.152676,	0.040424,	0.02195];
 end



%**************************************************************************
function LU =LU1996(cm,cf,cmf,V3,V2) %% es 1995
Bt13=(cm-cmf)/(cm+cmf);
Bt23=(cf-cmf)/(cf+cmf);

LU=cm*(1+2*(V2+V3)*(Bt23-Bt13*(V2+V3)/V3)/...
                   ((V2+V3)/V3 -Bt13*Bt23 - (V2+V3)*(Bt23-Bt13*(V2+V3)/V3)));
end

