clearvars
clear all 

ii=7; % 60
load Sk_Tk_30_mas_5_675_hasta_90_ENoP
delta1a=delta1; %2*ZZ(ii);
delta2a=delta2;
delta1=delta1a(ii);
delta2=delta2a(ii); 

Bi=[ 1.e12,	1,	0.1,	0.01,	0.001,	0.0001];
lam=[0.2,0.5,0.7,0.8,0.88,0.9,0.905,0.906899];
ro=[1001,101,11];
% ro1=[990.5, 99.498, 10.41953]; %%x1=cmf/cm;

cm=1;


Periodos2=(cos(theta(ii))+1i*sin(theta(ii)));
w1=1;
w2=Periodos2;
rm=min(sin(theta(ii)/2),0.5);
percolacion =pi*rm.^2./sin(theta(ii));

 %***********************************************************
% Metodo de Homogenizacion Asintotica
%*************************************************************
V = abs(w1)*abs(w2)*sin(theta(ii)); % *** Volumen de la celda periodica *******

%H1=(conj(delta1)*conj(w2)-conj(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
H1=(pi)/imag(conj(w1)*w2);
H2=((delta1)*conj(w2)-(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
S1=Sk(:,ii);
% *************************************************************************
for rr=1:length(ro) 
cf=ro(rr);
xas=cf/cm;
  for kk=1:length(Bi)
     for jj=1:length(lam)
     V2=lam(jj); % volumen de la fibra
     K=cf*Bi(kk); %% 10.41953;  %%%99.498; %%990.5;
     r1=sqrt(V*V2/pi);  % radio fibra 
   % *** Radio de la fibra más la mesophase ************************************
     X1=(xas - K*(xas - 1))/(xas + K*(xas + 1));
     J=eye(2)+X1*r1^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)];
      %CENo(jj,:,kk,rr) =coeficientes_efectivos(cm,V2,X1,J);
     CENo(jj,:,kk,rr) =ahm_imperfect_Paralelogramo_Cortas_2023_No6(cm,xas,K,V2,S1,J,r1);
     %CENo1(jj,:,kk,rr) = ahm_imperfect_Paralelogramo(cm,cf,K,V2,V,S1,H1,H2,35);
     CENo1(jj,:,kk,rr) = conductividad_efectiva_imperfect_orden_m(cm,xas,K,V2,S1,J,r1,42);
     CELu(jj,kk,rr)=LU1996(cm,cf,V2,K/cf);
     end
  end
end
% 
% %%% para comparar formula corta de Lu con O_0 nuestra
% Tabla=[];
% for kk=1:length(Bi)   
% Tabla=[Tabla,CENo(:,1,kk,1), CELu(:,kk,1)];
% %textoBi=[textoBi; num2str(Bi(kk))];
% end
%  disp(num2str([lam' Tabla],6))
% %  
%% %Tabla
 Lu95hex=Lu1995Tabla4;
disp('Comparacion Lu1995 con spring AHM')
% Configuración de formato
decimales_lam = 6;   % Decimales para lambda
decimales_val = 6;   % Decimales para valores AHM/Lu
ancho_col = 15;      % Ancho de cada columna
Bi = [1.e12, 1, 0.1, 0.01, 0.001, 0.0001]; % Valores de Bi
num_bi = length(Bi); % Debe ser par (6 en este caso)

for rr = 1:length(ro)
    fprintf('\rrho= %.4f\n', ro(rr))
    
    % Primera mitad (Bi1-Bi3)
    %fprintf('\nPrimera Parte (Bi1-Bi3):\n');
    % Primero: Imprimir los valores de Bi (centrados sobre AHM/Lu)
    disp('-----------------------------------------------------------------------------------------------------------')
    fprintf(['%', num2str(ancho_col), 's'], ' '); % Espacio para Lambda
    for kk = 1:num_bi/2
        fprintf(['%', num2str(ancho_col), '.1e%', num2str(ancho_col), '.1e'], Bi(kk), Bi(kk));
    end
    fprintf('\n');
    % Segundo: Imprimir encabezados (Lambda, AHM1, Lu1, AHM2, Lu2, ...)
    fprintf(['%', num2str(ancho_col), 's'], 'Lambda');
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
        fprintf(['%', num2str(ancho_col), '.1e%', num2str(ancho_col), '.1e'], Bi(kk), Bi(kk));
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
Bi = [1.e12, 1, 0.1, 0.01, 0.001, 0.0001]; % Valores de Bi
num_bi = length(Bi); % Debe ser par (6 en este caso)

for rr = 1:length(ro)
    fprintf('\rrho= %.4f\n', ro(rr))
     disp('-----------------------------------------------------------------------------------------------------------')
    % Primera mitad (Bi1-Bi3)
    %fprintf('\nPrimera Parte (Bi1-Bi3):\n');
    % Primero: Imprimir los valores de Bi (centrados sobre AHM/Lu)
    fprintf(['%', num2str(ancho_col), 's'], ' '); % Espacio para Lambda
    for kk = 1:num_bi/2
        fprintf(['%', num2str(ancho_col), '.1e%', num2str(ancho_col), '.1e'], Bi(kk), Bi(kk));
    end
    fprintf('\n');
    % Segundo: Imprimir encabezados (Lambda, AHM1, Lu1, AHM2, Lu2, ...)
    fprintf(['%', num2str(ancho_col), 's'], 'Lambda');
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
        fprintf(['%', num2str(ancho_col), '.1e%', num2str(ancho_col), '.1e'], Bi(kk), Bi(kk));
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
    for ii = 6:length(lam)
        fprintf(['%', num2str(ancho_col), '.', num2str(decimales_lam), 'f'], lam(ii));
        for kk = num_bi/2+1:num_bi
            fprintf(['%', num2str(ancho_col), '.', num2str(decimales_val), 'f'], CENo1(ii,1,kk,rr));
            fprintf(['%', num2str(ancho_col), '.', num2str(decimales_val), 'f'], Lu95hex(ii,kk,rr));
        end
        fprintf('\n');
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%
load S10_45_75_90_delta1_2; 
%%%%% primera celda 
ii=4; % 45
Bi=[0.:.00005:.1, .11:.01:.41, .5:.05:1]; 
lam=[.65,.64,.62, .6:-.1:.2];
lam1(1,:)=lam;
cm=1;
S1=Sk(:,ii);
Periodos2=(cos(theta(ii))+1i*sin(theta(ii)));
w1=1;
w2=Periodos2;
rm=min(sin(theta(ii)/2),0.5);
percolacion =pi*rm.^2./sin(theta(ii));
delta1=delta1a(ii);
delta2=delta2a(ii); 
 %***********************************************************
% Metodo de Homogenizacion Asintotica
%*************************************************************
V = abs(w1)*abs(w2)*sin(theta(ii)); % *** Volumen de la celda periodica *******
angulo(1,:)=['$\theta$ =' num2str(theta(ii)*180/pi,3),'$^o$'];
%H1=(conj(delta1)*conj(w2)-conj(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
H1=(pi)/imag(conj(w1)*w2);
H2=((delta1)*conj(w2)-(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
% *************************************************************************
for rr = 1:length(ro)
    cf=ro(rr);
    xas=cf/cm;
   for kk=1:length(Bi)
     for jj=1:length(lam)
     V2=lam(jj); % volumen de la fibra
     K=cf*Bi(kk); %% 10.41953;  %%%99.498; %%990.5;
    r1=sqrt(V*V2/pi);  % radio fibra 
   % *** Radio de la fibra más la mesophase ************************************

     X1=(xas - K*(xas - 1))/(xas + K*(xas + 1));
     J=eye(2)+X1*r1^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)];
 
    CENo(jj,:,kk,rr,1) =ahm_imperfect_Paralelogramo_Cortas_2023_No6(cm,xas,K,V2,S1,J,r1);
 %CENo(jj,:,kk) = ahm_imperfect_Paralelogramo(cm,cf,K,V2,V,S1,H1,H2,15);
     end
   end
end
 %%%%% Segunda  celda 
ii=11; % 75
%lam=[.81,.8:-.1:.2];
lam=[.75, .72 .7:-.1:.2];
lam1(2,:)=lam;
cm=1;
xas=cf/cm;
S1=Sk(:,ii);
Periodos2=(cos(theta(ii))+1i*sin(theta(ii)));
w1=1;
w2=Periodos2;
rm=min(sin(theta(ii)/2),0.5);
percolacion =pi*rm.^2./sin(theta(ii));
delta1=delta1a(ii);
delta2=delta2a(ii); 
 %***********************************************************
% Metodo de Homogenizacion Asintotica
%*************************************************************
V = abs(w1)*abs(w2)*sin(theta(ii)); % *** Volumen de la celda periodica *******
angulo(2,:)=['$\theta$ =' num2str(theta(ii)*180/pi,3),'$^o$'];
%H1=(conj(delta1)*conj(w2)-conj(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
H1=(pi)/imag(conj(w1)*w2);
H2=((delta1)*conj(w2)-(delta2)*conj(w1))/(w1*conj(w2)-w2*conj(w1));
h11=real(H1);
h21=imag(H1);
h12=real(H2);
h22=imag(H2);
% *************************************************************************
for rr = 1:length(ro)
    cf=ro(rr);
    xas=cf/cm;
for kk=1:length(Bi)
for jj=1:length(lam)
      V2=lam(jj); % volumen de la fibra
      K=cf*Bi(kk); %% 10.41953;  %%%99.498; %%990.5;
      r1=sqrt(V*V2/pi);  % radio fibra 
      % *** Radio de la fibra más la mesophase ************************************

      X1=(xas - K*(xas - 1))/(xas + K*(xas + 1));
      J=eye(2)+X1*r1^2*[(h11+h12),(h21-h22);(-h21-h22),(h11-h12)];
 
     CENo(jj,:,kk,rr,2) =ahm_imperfect_Paralelogramo_Cortas_2023_No6(cm,xas,K,V2,S1,J,r1);
     %CENo(jj,:,kk) = ahm_imperfect_Paralelogramo(cm,cf,K,V2,V,S1,H1,H2,15);
     end
   end
end

matriz=[ 'cm= ',num2str(cm,5)];
fibra=[ 'cf= ',num2str(cf,5)];
rho=['$\rho $','=',num2str(cf/cm,5)];

Modulo=[ '|w2|= ',num2str(abs(w2),5)];
volmefiber=[ '$V_2= $',num2str(V2,3)];
% disp(Modulo)
% disp(angulo)
% disp(matriz)
% disp(fibra)
% disp(mesofase)

 % Crear el texto dinámico para la leyenda
% texto_leyenda = sprintf('Vector: [%s]', num2str(Bi));
ejey = ['$\hat{\kappa}_{11}/\kappa_1$'; '$\hat{\kappa}_{22}/\kappa_1$'; '$\hat{\kappa}_{12}/\kappa_1$'];
% lam = V2masV3;

% Definir una paleta de colores (ej. usando 'lines', 'parula', 'jet', o personalizada)
colors = lines(length(Bi));  % 'lines' genera colores distintos, también puedes usar 'jet', 'hsv', etc.
colores=[[0,0,1];[.6,.6,0];[0,1,0];[1,0,0];[0,1,1];[1,0,1];[0,0,0];[0.5,0.5,0.5];[0.5,0.5,1];[0.8,0.8,0];[0.1,0.8,0];...
         [.8,0,0];[0,.8,.8];[.1,.1,0];[0,.1,0];[.1,0,0]; [0.9, 0.7, 0.5 ]; [0.75, 0.5, 0.5];[0.25, 0.25, 0.5];...
          [0.75, 0.75, 0];[0.5, 0, 0.5]];
 incisos=['a)'; 'b)'; 'c)'; 'd)'; 'e)';'f)'];
for rr=3:length(ro)
figure(rr)
for ii=1:2
     lam=lam1(ii,:);
%set(gcf, 'Position', [100 100 1200 400]); % Ajusta el tamaño de la figura   
   for jj = 1:3
    subplot(2,6, 3*(ii-1)+jj)
    cm=1; cf=ro(rr);
    rho1=[ '$\rho$','=',num2str(cf/cm,5)];
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
    h(kk) = plot(Bi, CENo_z, ...
                'LineStyle', linestyle, ...
                'LineWidth', linewidth, ...
                'Color', colores(kk,:));
    legendEntries{kk} = ['$V_f = $' num2str(lam(kk),3)];
    end
 
    hold off
    
    % Etiquetas de ejes
     ylabel(ejey(jj,:), 'fontsize', 14, 'interpreter', 'latex', 'FontWeight', 'bold');
    xlabel(['$Bi$'], 'fontsize', 10, 'interpreter', 'latex', 'FontWeight', 'bold');
    
    if jj==1
    % Leyenda
    legend(h, legendEntries, 'FontSize', 10, 'Location', 'best', 'interpreter', 'latex', 'FontWeight', 'bold','NumColumns',4);
    title( [rho1,', ', angulo],'fontsize', 12, 'interpreter', 'latex')
    end
    
    % Cuadrícula y marco
    grid on
    box on
    set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.7)
    
    % Anotaciones de texto
    text_pos_x = min(Bi) + 0.1*(max(Bi)-min(Bi));
    
    % Textos fijos
    %title( [rho1,', ',rho2,', ',  angulo],'fontsize', 10, 'interpreter', 'latex')
    %title( [angulo],'fontsize', 12, 'interpreter', 'latex', 'FontWeight', 'bold');
%     text(text_pos_x, CENo(1,jj,2,ii), rho1,'fontsize', 14, 'interpreter', 'latex')
%     text(text_pos_x, CENo(1,jj,5,ii), rho2, 'fontsize', 14, 'interpreter', 'latex')
%     text(text_pos_x, CENo(1,jj,7,ii), angulo, 'fontsize', 14, 'interpreter', 'latex')
%     text(.08, max(ylim), incisos(jj,:), 'fontsize', 12, 'interpreter', 'latex')
    % Ajuste de ejes
    xlim([min(Bi) max(Bi)])
    %ylim auto
    
    %%%%%%%%%%%%%
    end
   end
end

for rr=3:length(ro)
figure(rr)
for ii=1:2
     lam=lam1(ii,:);
%set(gcf, 'Position', [100 100 1200 400]); % Ajusta el tamaño de la figura   
   for jj = 1:3
    subplot(2,6, 3*(ii-1) + jj + 6)
    cm=1; cf=ro(rr);
     rho1=[ '$\rho$','=',num2str(cf/cm,5)];
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
    h(kk) = plot(Bi, CENo_z, ...
                'LineStyle', linestyle, ...
                'LineWidth', linewidth, ...
                'Color', colores(kk,:));
    legendEntries{kk} = ['$V_f = $' num2str(lam(kk),3)];
    end
    hold off
    
    % Etiquetas de ejes
     ylabel(ejey(jj,:), 'fontsize', 14, 'interpreter', 'latex', 'FontWeight', 'bold');
   xlabel(['$Bi$'], 'fontsize', 10, 'interpreter', 'latex', 'FontWeight', 'bold');
    
    
    if jj==1
    % Leyenda
    legend(h, legendEntries, 'FontSize', 10, 'Location', 'best', 'interpreter', 'latex', 'FontWeight', 'bold','NumColumns',4);
    title( [rho1,', ',  angulo],'fontsize', 12, 'interpreter', 'latex')
    end
    
    % Cuadrícula y marco
    grid on
    box on
    set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.7)
    
    % Anotaciones de texto
    text_pos_x = min(Bi) + 0.1*(max(Bi)-min(Bi));
    
    % Textos fijos
    %title( [rho1,', ',rho2,', ',  angulo],'fontsize', 10, 'interpreter', 'latex')
    %title( [ angulo],'fontsize', 12, 'interpreter', 'latex', 'FontWeight', 'bold');
%     text(text_pos_x, CENo(1,jj,2,ii), rho1,'fontsize', 14, 'interpreter', 'latex')
%     text(text_pos_x, CENo(1,jj,5,ii), rho2, 'fontsize', 14, 'interpreter', 'latex')
%     text(text_pos_x, CENo(1,jj,7,ii), angulo, 'fontsize', 14, 'interpreter', 'latex')
%     text(.08, max(ylim), incisos(jj,:), 'fontsize', 12, 'interpreter', 'latex')
    % Ajuste de ejes
    xlim([0 .2])
    %ylim auto
    end
   end
end

% Bi=[0.0001,  0.01,	0.1,	0.2]; 
figure(4)
for ii=1:2
     lam=lam1(ii,:);
for rr=1:length(ro)
    
%set(gcf, 'Position', [100 100 1200 400]); % Ajusta el tamaño de la figura   

    subplot(2,6, 3*(ii-1) + rr)
     cm=1; cf=ro(rr);
    rho1=[ '$\rho$','=',num2str(cf/cm,5)];
    angulo=['$\theta$ =' num2str(theta(4*ii+1)*180/pi,3),'$^o$'];
    hold on
    % Prealocación manual para handles y leyenda
    h = zeros(length(lam), 1);

    legendEntries = cell(length(lam), 1);
    hold on
    
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
    
    % Extracción compatible con versiones antiguas
        CENo_z1 = reshape(CENo(kk, 1, :,rr,ii), [], 1); % Equivalente a squeeze
        CENo_z2 = reshape(CENo(kk, 2, :,rr,ii), [], 1); % Equivalente a squeeze
        % Graficar con estilo
        h(kk) = plot(Bi, CENo_z1./CENo_z2 , '-', ...
                'LineStyle', linestyle, ...
                'LineWidth', linewidth, ...
                'Color', colores(kk,:));
    legendEntries{kk} = ['$V_f = $' num2str(lam(kk),3)];
    end
    hold off
  
    
    % Etiquetas de ejes
     ylabel('$\hat{\kappa}_{11}/\hat{\kappa}_{22}$', 'fontsize', 14, 'interpreter', 'latex', 'FontWeight', 'bold');
   xlabel(['$Bi$'], 'fontsize', 10, 'interpreter', 'latex', 'FontWeight', 'bold');
    
    
    if rr==1
    % Leyenda
    legend(h, legendEntries, 'FontSize', 10, 'Location', 'best', 'interpreter', 'latex', 'FontWeight', 'bold','NumColumns',4);
    end
    
    title( [rho1,', ',  angulo],'fontsize', 12, 'interpreter', 'latex')
    
    
    % Cuadrícula y marco
    grid on
    box on
    set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.7)
    
    % Anotaciones de texto
    text_pos_x = min(Bi) + 0.1*(max(Bi)-min(Bi));
    
    % Ajuste de ejes
    xlim([1/(cf-1)+.005,  max(Bi)])
  
    %ylim auto
   
    % Textos fijos
%     text(text_pos_x, 0.9, rho,'fontsize', 14, 'interpreter', 'latex')
%     text(text_pos_x, 0.7, angulo, 'fontsize', 14, 'interpreter', 'latex')
    
      subplot(2,6, 3*(ii-1) + rr+6)
     cm=1; cf=ro(rr);
    rho1=[ '$\rho$','=',num2str(cf/cm,5)];
    angulo=['$\theta$ =' num2str(theta(4*ii+1)*180/pi,3),'$^o$'];
    hold on
    % Prealocación manual para handles y leyenda
    h = zeros(length(lam), 1);

    legendEntries = cell(length(lam), 1);
    hold on
    
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
    
    % Extracción compatible con versiones antiguas
        CENo_z1 = reshape(CENo(kk, 1, :,rr,ii), [], 1); % Equivalente a squeeze
        CENo_z2 = reshape(CENo(kk, 2, :,rr,ii), [], 1); % Equivalente a squeeze
        % Graficar con estilo
        h(kk) = plot(Bi, CENo_z1./CENo_z2 , '-', ...
                'LineStyle', linestyle, ...
                'LineWidth', linewidth, ...
                'Color', colores(kk,:));
    legendEntries{kk} = ['$V_f = $' num2str(lam(kk),3)];
    end
    hold off
  
    
    % Etiquetas de ejes
     ylabel('$\hat{\kappa}_{11}/\hat{\kappa}_{22}$', 'fontsize', 14, 'interpreter', 'latex', 'FontWeight', 'bold');
   xlabel(['$Bi$'], 'fontsize', 10, 'interpreter', 'latex', 'FontWeight', 'bold');
    
    
    if rr==1
    % Leyenda
    legend(h, legendEntries, 'FontSize', 10, 'Location', 'best', 'interpreter', 'latex', 'FontWeight', 'bold','NumColumns',4);
    end
    
    title( [rho1,', ',  angulo],'fontsize', 12, 'interpreter', 'latex')
    
    
    % Cuadrícula y marco
    grid on
    box on
    set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.7)
    
    % Anotaciones de texto
    text_pos_x = min(Bi) + 0.1*(max(Bi)-min(Bi));
    
    % Ajuste de ejes
    
      xlim([0, 2/(cf-1)])
    %ylim auto
   
    % Textos fijos
%     text(text_pos_x, 0.9, rho,'fontsize', 14, 'interpreter', 'latex')
%     text(text_pos_x, 0.7, angulo, 'fontsize', 14, 'interpreter', 'latex')
    
   
  end
end




%%%%%%%%%%%%%


 
%%% FUNCIONES AUXILIARES

%%% orden del sistema No=1
function RR1 =ahm_imperfect_Paralelogramo_Cortas_2023_No1(cm,xas,K,V2,S1,J,r1)
              
S4=S1(4);S6=S1(6);S8=S1(8);S10=S1(10);S12=S1(12);S14=S1(14);S16=S1(16);S18=S1(18);
S20=S1(20);S22=S1(22);S24=S1(24);

X1=(xas - K*(xas - 1))/(xas + K*(xas + 1));
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
function RR2 =ahm_imperfect_Paralelogramo_Cortas_2023_No2(cm,xas,K,V2,S1,J,r1)

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
function RR3 =ahm_imperfect_Paralelogramo_Cortas_2023_No3(cm,xas,K,V2,S1,J,r1)

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
function RR4 =ahm_imperfect_Paralelogramo_Cortas_2023_No4(cm,xas,K,V2,S1,J,r1)

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
W21=7*sqrt(15)*X5*r1^8*L8;     W22=126*X5*r1^10*L10;           W23=66*sqrt(35)*X5*r1^12*L12;    W24=143*sqrt(45)*X5*r1^14*L14;  
W31=12*sqrt(21)*X7*r1^10*L10;  W32=66*sqrt(35)*X7*r1^12*L12;   W33=1716*X7*r1^14*L14;           W34=715*sqrt(63)*X7*r1^16*L16;    
W41=55*sqrt(3)*X9*r1^12*L12;  W42=143*sqrt(45)*X9*r1^14*L14;   W43=715*sqrt(63)*X9*r1^16*L16;   W44=24310*X9*r1^18*L18;           


W=[W11, W12, W13, W14; W21, W22, W23, W24;W31, W32, W33, W34; W41, W42, W43, W44];

AP2=eye(8)+W;
PP=-X1*N1*inv(AP2)*N2;  
JK=J+PP; % Matriz Z
% disp('corto');
Z=JK;

RR4 =coeficientes_efectivos(cm,V2,X1,Z);

end


%%% orden del sistema No=5
function RR5 =ahm_imperfect_Paralelogramo_Cortas_2023_No5(cm,xas,K,V2,S1,J,r1)

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
W21=7*sqrt(15)*X5*r1^8*L8;     W22=126*X5*r1^10*L10;           W23=66*sqrt(35)*X5*r1^12*L12;    W24=143*sqrt(45)*X5*r1^14*L14;     W25=273*sqrt(55)*X5*r1^16*L16;
W31=12*sqrt(21)*X7*r1^10*L10;  W32=66*sqrt(35)*X7*r1^12*L12;   W33=1716*X7*r1^14*L14;           W34=715*sqrt(63)*X7*r1^16*L16;    W35=1768*sqrt(77)*X7*r1^18*L18; 
W41=55*sqrt(3)*X9*r1^12*L12;  W42=143*sqrt(45)*X9*r1^14*L14;   W43=715*sqrt(63)*X9*r1^16*L16;   W44=24310*X9*r1^18*L18;           W45=8398*sqrt(99)*X9*r1^20*L20; 
W51=26*sqrt(33)*X11*r1^14*L14; W52=273*sqrt(55)*X11*r1^16*L16; W53=1768*sqrt(77)*X11*r1^18*L18; W54=8398*sqrt(99)*X11*r1^20*L20; W55=352716*X11*r1^22*L22;

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
function RR6 =ahm_imperfect_Paralelogramo_Cortas_2023_No6(cm,xas,K,V2,S1,J,r1)

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
W21=7*sqrt(15)*X5*r1^8*L8;     W22=126*X5*r1^10*L10;           W23=66*sqrt(35)*X5*r1^12*L12;    W24=143*sqrt(45)*X5*r1^14*L14;     W25=273*sqrt(55)*X5*r1^16*L16;
W31=12*sqrt(21)*X7*r1^10*L10;  W32=66*sqrt(35)*X7*r1^12*L12;   W33=1716*X7*r1^14*L14;           W34=715*sqrt(63)*X7*r1^16*L16;    W35=1768*sqrt(77)*X7*r1^18*L18; 
W41=55*sqrt(3)*X9*r1^12*L12;  W42=143*sqrt(45)*X9*r1^14*L14;   W43=715*sqrt(63)*X9*r1^16*L16;   W44=24310*X9*r1^18*L18;           W45=8398*sqrt(99)*X9*r1^20*L20; 
W51=26*sqrt(33)*X11*r1^14*L14; W52=273*sqrt(55)*X11*r1^16*L16; W53=1768*sqrt(77)*X11*r1^18*L18; W54=8398*sqrt(99)*X11*r1^20*L20; W55=352716*X11*r1^22*L22;
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


function [kappa_eff] = conductividad_efectiva_imperfect_orden_m(kappa1,rho,K,V2,Sk,J1,R1,m)
                    
%     rho = kappa2 / kappa1;
    
    % Calculate chi values up to required order
    %%X9=(9*xas - K*(xas - 1))/(9*xas + K*(xas + 1));
    max_p = 6*(m) + 1; % We need up to chi_{2(m)+1}
    chi = zeros(1, max_p);
    for p = 1:2:max_p
        numerator = p*rho-K*(rho-1);
        denominator = p*rho+K*(rho+1);
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
    kappa_eff = coeficientes_efectivos(kappa1,V2,chi(1),Z);
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
for n=1:2:71
   for k=1:2:71

Wnk(n,k)=(S(n+k))*(r^(n+k))*((factorial(n+k-1)/(factorial(n-1)*factorial(k-1)))/...
               sqrt(n*k));

   end
end
RWnk=real(Wnk);
IWnk=imag(Wnk);
di=eye(70); 

% Para hallar los Bp

for s=1:2:90,

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

 function Lu95hex=Lu1995Tabla4
 
 Lu95hex(:,:,1)=[ 1.49876        1.49751        1.48649        1.39046              1       0.718793       0.666665;
                  2.99671        2.98873        2.91896        2.38158              1       0.419058       0.332809;
                  5.77799        5.74463         5.4618        3.69406              1       0.268559       0.172083;
                  9.85697        9.75231        8.90664         4.8678              1       0.201035       0.100415;
                  24.7331        23.9164        18.6074        6.41883              1       0.145317       0.039292;
                  52.4503        47.5269        27.6829        6.96727              1       0.127947       0.017784;
                  96.96        75.6144        32.2877        7.11943              1        0.12199        0.00887;
                   0          113.2         34.648          7.179              1          0.119              0];
     
          
  Lu95hex(:,:,2)=[ 1.48781        1.47596        1.38298              1       0.719188       0.672238       0.666665;
                  2.92736         2.8542        2.34441              1       0.419732       0.341773       0.332809;
                   5.4991        5.21554        3.59391              1       0.269325       0.181984       0.172083;
                   9.03672        8.24214        4.68814              1       0.201837       0.110699       0.100415;
                   19.787        15.7537        6.10023              1       0.146162       0.050142       0.039292;
                   32.9927        21.5627        6.59024              1        0.12882       0.029139       0.017784;
                   43.95        24.1169        6.72533              1        0.12288         0.0206        0.00887;
                   55        25.3338        6.77812              1           0.12              0              0];
                  
    Lu95hex(:,:,3)=[1.4        1.32143              1       0.723076       0.672726          0.667276         0.666665;
                    2.43058        2.05944              1       0.426391        0.34256       0.333788       0.332809;
                    3.84248        2.89114              1       0.276891       0.182855       0.173164       0.172083;
                    5.18184         3.5192              1       0.209772       0.111604       0.101537       0.100415;
                    7.1887        4.21427              1       0.154523       0.051097       0.040475       0.039292;
                    8.05223         4.4291              1       0.137454       0.030137       0.019024       0.017784;
                    8.3264        4.48615              1        0.13169        0.02166        0.01015        0.00887;
                    8.44102         4.5082              1          0.129              0              0              0];
 end
 
 
 
%**************************************************************************
function LU =LU1996(cm,cf,V2,K)
xas=cf/cm;
LU=cm*(1+ 2*V2*(K*(xas-1)-1)/...
      (K*(xas+1)+1- V2*(K*(xas-1)-1)));
end


