% Definir la señal
N=256;
n=0:N-1;
xn =  3*cos(2*pi*n/10);

% Encontrar la potencia de 2 más cercana y rellenar con ceros si es necesario
N = length(xn); %Longitud de la señal de entrada    
N_potencia_2 = 2^nextpow2(N); %Potencia de 2 más cercana al N de la señal de entrada
xn_new = [xn, zeros(1, N_potencia_2 - N)]; %Creamos una nueva señal seguida de 0 de ser necesarios para completar la potencia más cercana
%tic toc para la función propia
tic
senal_t = FFTButterfly(xn_new);
tiempoFuncion = toc;
%Tic toc para función propia de MATLAB
tic
fft(xn_new)
tiempoFFTMatlab = toc;
%Tic toc para DFT
tic
xk= zeros(1,N);
for k=0:N-1
    for n=0:N-1
    xk(k+1) = xk(k+1) + xn(n+1) * exp(-1i * 2*pi * k * n / N);
    end
end
tiempoDFT = toc;
fprintf('Resultados:\n');
disp(senal_t.');
%% Comparación con DFT del Deber 4
disp(['Tiempo transcurrido para la función realizada: ', num2str(tiempoFuncion), ' segundos']);
disp(['Tiempo transcurrido para la función DFT del Deber 4: ', num2str(tiempoDFT), ' segundos']);
%% Comparación de FFT MATLAB vs Propia
disp(['Tiempo transcurrido para la función realizada: ', num2str(tiempoFuncion), ' segundos']);
disp(['Tiempo transcurrido para la función propia de MATLAB: ', num2str(tiempoFFTMatlab), ' segundos']);
%% Función
function a = FFTButterfly(x)
    N = length(x); % Longitud de la señal
    M = log2(N); % Número total de etapas
    % Reordenamiento de la señal
    x = bitrevorder(x);
    for etapa = 1:M % Iterar sobre cada etapa
        paso = 2^(etapa-1); % Tamaño del paso, se duplica en cada etapa
        salto = 2^etapa; % Salto de índices al iterar sobre bloques de la señal en cada etapa
        w_m = exp(-2i * pi / salto); % Factor de giro para esta etapa

        for inicio = 1:salto:N % Iterar sobre cada bloque en esta etapa, el salto representa el tamaño del bloque
            w = 1; % Inicializar el factor de giro
            for k = inicio:inicio+paso-1 % Iterar sobre cada par de puntos en el bloque actual
                i_par = k; % Índice par
                i_impar = k + paso; % Índice impar

                % Operaciones de mariposa
                par = x(i_par); %Guardamos el valor del punto par
                impar = x(i_impar) * w; %Se multiplica el valor del punto impar por un factor de giro w, que se calcula para esta etapa específica
                x(i_par) = par + impar; %Actualizar los valores de los puntos par e impar con el resultado de las operaciones de mariposa
                x(i_impar) = par - impar;

                % Actualizar el factor de giro
                w = w * w_m;
            end
        end
    end
    a = x;
end