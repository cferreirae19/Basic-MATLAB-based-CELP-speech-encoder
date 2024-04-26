
function [SNRseg,SNRm,m] = SNRS(x,xq,L)

% -Funcionalidad:
%    ·Calcular la relación señal a ruido por segmentos de una señal dada
% -Parámetros de entrada:
%    ·x: Señal original
%    ·xq: Señal cuantificada
%    ·L: Longitud del segmento
% -Parámetros de salida:
%    ·SNRm: Relación señal a ruido por tramas (vector)
%    ·SNRseg: relación señal a ruido segmental
%    ·m: Vector de referencia en tiempo discreto para S

% Dividir señal original en tramas de longitud L

n = numel(x);           % 'numel' devuelve el nº de elementos en un array
xn = mat2cell(x,diff([0:L:n-1,n]));            % 'mat2cell' convierte una matriz en una matriz de celdas (estas celdas contienen submatrices)

% Dividir señal cuantificada en tramas de longitud L

xn_q = mat2cell(xq,diff([0:L:n-1,n]));

% Cálculo de la SNR por tramas

c = numel(xn);          % Número de celdas del array de celdas xn (es el mismo número que para xn_q)
k = 0:c-1;
for i = 1:c
    SNRm(i) = SNR(xn{i,1},xn_q{i,1});           % Cálculo de SNRm(i), que es cada una de las tramas de la SNR
end

SNRseg = mean(SNRm);            %Cálculo de la relación señal a ruido segmental (valor medio de la SNR por tramas)

m = L:L:n;          % Índices de x correspondientes a la última muestra de cada trama de la SNRm
if mod(n,L) ~= 0            %  Si la última trama es de menor longitud que L
    m(end+1) = n;           % Añadimos un último índice (el correspondiente al último elemento de la señal)           
end