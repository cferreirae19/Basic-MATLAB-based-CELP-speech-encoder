function [sh,B,G,AK,Tv,indv,bits_muestra,ep,ead,eest,zaux1,zaux2]=celp_basico(s,Ltrama,Lsubtrama,p)

% sh, señal de voz reconstruida
% B y G vectores con las ganancias adaptativa y estocástica de las
% distintas subtramas
% Ak, matriz con los LPC de cada trama
% bits_muestra, número de bits por muestra.
% Tv, vector con los distintos valores de T en la excitación adaptativa.
% indv, vector con los distintos índices seleccionados de la excitación estocástica.
%ep error de predicción: ead parte adaptativa, eest parte estocástica

% s, señal de voz
% Ltrama y Lsubtrama, longitudes de trama y subtrama.
% p, orden de predicción.

s_tr = s.';

sh = [];
ep = [];
ead = [];
eest = [];
B = [];
G = [];
AK = [];
Tv = [];
indv = [];

window = rectwin(Ltrama);
wshift = Ltrama;
zaux1 = [];
zaux2 = [];
adapt = zeros(1,3*Lsubtrama);

% Inicialización de la biblioteca estocástica
M = 512;
N = Lsubtrama;
rng(1);
v = randn(M,N);

for k=1:length(s_tr)/Ltrama

    [AK(k,:),~,~,~,~] = speech2lpc(s(((k-1)*Ltrama+1):k*Ltrama),p,window,wshift);   % Estimamos los coeficientes del filtro A(z) para la subtrama actual

    for i=1:4   % Para cada una de las 4 subtramas que conforman 1 trama

        max_index_t = 0;
        max_gain_b = 0;
        max_index_l = 0;
        max_gain_g = 0;
        u0_def = zeros(1,Lsubtrama);
        uf_def = zeros(1,Lsubtrama);
        min_Et = 1000;   % Inicializamos las energías de error a un valor bastante alto para que en alguna iteración se pueda conseguir un valor menor
        min_El = 1000;
        y1l_ret = zeros(1,Lsubtrama);
        y2t_ret = zeros(1,Lsubtrama);

        %=========================================%
        %==========BIBLIOTECA ADAPTATIVA==========%
        %=========================================%

        % Probaremos con todos los valores de t los siguientes pasos:
            % 1: Extraer candidato de la biblioteca adaptativa (d20t)
            % 2: Filtrar d20t con un filtro 1/A(z) para obtener y20t
            % 3: Calcular la ganancia bt
            % 4: Obtener la energía del error
            % 5: Si es mayor o igual que el mínimo almacenado, repetir 
            %    desde el paso 1 con el siguiente valor de t (t++)
            % 6: Si es menor, guardar el valor actual del índice (t), la
            %    ganancia (bt) y la energía (u0) y repetir desde el paso 1
            %    con el siguiente valor de t (t++)

        for t=1:(2*Lsubtrama+1)

            d20t = adapt(t:t+Lsubtrama-1);  % Extraemos un candidato de la biblioteca adaptativa
            [y20t] = filter(1,AK(k,:),d20t,zaux1);    % Filtramos el candidato para obtener y20t
            y20t_tr = y20t.';
            bt = (s_tr(((k-1)*Ltrama+(i-1)*Lsubtrama+1):((k-1)*Ltrama+i*Lsubtrama))*y20t_tr)/(y20t*y20t_tr+eps);   % Cálculo de la ganancia
            y2t = bt*y20t;  % Multiplicamos por la ganancia (en ese orden para evitar errores con los tamaños de las matrices al hacer el producto)
            u0 = s_tr(((k-1)*Ltrama+(i-1)*Lsubtrama+1):((k-1)*Ltrama+i*Lsubtrama)) - y2t;
            u0_tr = u0.';
            Et = u0*u0_tr;  % Cálculo de la energía del error

            if (Et<min_Et)   % Si la energía del error es más pequeña que la menor encontrada hasta ahora
                                
                max_index_t = t;  % Guardamos la ganancia y el índice de la iteración actual porque son los valores que minimizan la energía del error
                max_gain_b = bt;
                min_Et = Et;    % Guardamos la energía del error para poder compararla con la de futuras iteraciones y saber cual es menor
                u0_def = u0;
                y2t_ret = y2t;

            end

        end
        
        d20t = adapt(max_index_t:max_index_t+Lsubtrama-1);  % Extraemos el candidato de la biblioteca adaptativa
        [y20t,zaux1] = filter(1,AK(k,:),d20t,zaux1);    % Filtramos el candidato para obtener y20t
        d2t=d20t*max_gain_b; %Obtener error adaptativo

        % Añadimos el valor de la ganancia bt al vector de retorno
        B = [B,max_gain_b];

        % Calculamos el valor del parámetro T y lo añadimos al vector de retorno
        T = 3*Lsubtrama - max_index_t + 1;
        Tv = [Tv,T];

        %==========================================%
        %==========BIBLIOTECA ESTOCÁSTICA==========%
        %==========================================%

        % Probaremos con todos los valores de l los siguientes pasos:
            % 1: Extraer candidato de la biblioteca estocástica (vl)
            % 2: Filtrar vl con un filtro 1/A(z) para obtener y10t
            % 3: Calcular la ganancia gl
            % 4: Obtener la energía del error
            % 5: Si es mayor o igual que el mínimo almacenado, repetir
            %    desde el paso 1 con el siguiente valor de l (l++)
            % 6: Si es menor, guardar el valor actual del índice (l), la
            %    ganancia (gl) y la energía (uf) y repetir desde el paso 1
            %    con el siguiente valor de l (l++)

        for l=1:length(v)

            vl = v(l,:);  % Extraemos un candidato de la biblioteca estocástica
            [y10l] = filter(1,AK(k,:),vl,zaux2);   % Filtramos el candidato para obtener y10l
            y10l_tr = y10l.';
            gl = (u0_def*y10l_tr)/(y10l*y10l_tr+eps);    % Cálculo de la ganancia
            y1l = gl*y10l;    % Multiplicamos por la ganancia (en ese orden para evitar errores con los tamaños de las matrices al hacer el producto)
            uf = u0_def - y1l;
            uf_tr = uf.';
            El = uf*uf_tr;      % Cálculo de la energía del error

            if (El<min_El)  % Si la energía del error es más pequeña que la menor encontrada hasta ahora

                max_index_l = l;  % Guardamos la ganancia y el índice de la iteración actual porque son los valores que minimizan la energía del error
                max_gain_g = gl;
                min_El = El;    % Guardamos la energía del error para poder compararla con la de futuras iteraciones y saber cual es menor
                uf_def = uf;
                y1l_ret = y1l;

            end
        
        end

        vl = v(max_index_l,:);  % Extraemos el candidato de la biblioteca estocástica
        [y10l,zaux2] = filter(1,AK(k,:),vl,zaux2);   % Filtramos el candidato para obtener y10l
        d1l=vl*max_gain_g; %Obtener error estocástico

        % Añadimos el valor de la ganancia gl al vector de retorno
        G = [G,max_gain_g];

        % Añadimos el valor del índice l al vector de retorno
        indv = [indv,max_index_l];

        % Antes de acabar con la subtrama, actualizamos el valor de d y lo añadimos a la biblioteca adaptativa
        d = d1l + d2t;
        adapt_1_new = adapt((Lsubtrama+1):2*Lsubtrama);
        adapt_2_new = adapt((2*Lsubtrama+1):3*Lsubtrama);
        adapt = [adapt_1_new,adapt_2_new,d];

        % Calculamos la subtrama correspondiente para el error y la señal reconstruida.
        enew = d2t + d1l;
        ep = [ep,enew];
        ead = [ead, d2t];
        eest = [eest, d1l];
        
        sh_new = y1l_ret + y2t_ret;
        sh = [sh,sh_new];

    end

end

% Cálculo del número de bits por muestra
bits_B = 8*length(B);
bits_G = 8*length(G);
bits_indv = 9*length(indv);
bits_Tv = 7*length(Tv);
bits_AK = 3*size(AK,1);
bits_muestra = (bits_B+bits_G+bits_indv+bits_Tv+bits_AK)/length(s);

end