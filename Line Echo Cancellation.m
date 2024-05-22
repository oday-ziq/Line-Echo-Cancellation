load('path.mat'); 

figure;
plot(path);  
xlabel('Sample');
ylabel('Amplitude');
title('Impulse Response of Echo Path');

figure;
fs = 8000;  
N = length(path);
f = linspace(0, fs, N);
mag = abs(fft(path));
mag_dB = 20*log10(mag);
plot(f(1:N/2), mag_dB(1:N/2));  
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Frequency Response of Echo Path');
%c
load('css.mat'); 
N1 = 5600*5;
figure;
plot(css);
xlabel('Sample');
ylabel('Amplitude');
title('Composite Source Signal (CSS)');
psd = css .* css;
figure;
plot(abs(fft(psd)))
figure;
repeated_signal = repmat(css, 1, 5);
power = repeated_signal .* repeated_signal;
spower = abs(fft(power));
plot(repeated_signal)
figure;
y = conv(repeated_signal,path);
y1 = y .* y;
plot(y)
pi = 10*log10((1/N1)*sum(power));
po = 10*log10((1/N1)*sum(y1));
ERL = pi - po;

path_mat = load('path.mat');
[h, w] = freqz(path_mat.path(1,:), 1, 8000);
css_mat = load('css.mat');
[Pxx, f] = periodogram(css_mat.css(1,:), [], [], 8000);

far_end = repmat(css_mat.css(1,:), 1, 10);
echo = conv(far_end, path_mat.path(1,:));
w_0 = zeros(1, 128);
echo = echo(1:length(far_end));
N = length(far_end);
M = 128;
mu = 0.25;
w = w_0;
e = zeros(1, N);
eps = 1e-6;

for n = M:N
    x_n = fliplr(far_end(n-M+1:n)); 
    y_n = dot(w, x_n); 
    e(n) = echo(n) - y_n;
    w = w + mu * e(n) * x_n / (dot(x_n, x_n) + eps); 
end

figure;
plot(far_end);
title('Far-End Signal');

figure;
plot(echo);
title('Echo Signal');

figure;
plot(e);
title('Error Signal');

[h_est, w_est] = freqz(w, 1, 8000);

[h_path, w_path] = freqz(path_mat.path(1,:), 1, 8000);

figure;
plot(w_path, 20*log10(abs(h_path)), 'b', 'LineWidth', 1.5);
hold on;
plot(w_est, 20*log10(abs(h_est)), 'r', 'LineWidth', 1.5);
title('Amplitude Response');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('Given FIR System', 'Estimated FIR Channel');

figure;
plot(w_path, unwrap(angle(h_path)), 'b', 'LineWidth', 1.5);
hold on;
plot(w_est, unwrap(angle(h_est)), 'r', 'LineWidth', 1.5);
title('Phase Response');
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
legend('Given FIR System', 'Estimated FIR Channel');

w = w_0;
e = zeros(1, N);
eps = 1e-6;

for n = M:N
    x_n = fliplr(far_end(n-M+1:n));
    y_n = dot(w, x_n); 
    e(n) = echo(n) - y_n; 
    w = w + mu * e(n) * x_n / (dot(x_n, x_n) + eps);
end


figure;
plot(path_mat.path(1,:), 'b', 'LineWidth', 1.5);
hold on;
plot(w, 'r', 'LineWidth', 1.5);
title('Impulse Response and its estimate');
xlabel('Sample');
ylabel('Amplitude');
legend('Impulse Response', 'Estimated Impulse Response');


path_mat = load('path.mat');
css_mat = load('css.mat');
far_end = repmat(css_mat.css(1,:), 1, 10);
echo = conv(far_end, path_mat.path(1,:));
echo = echo(1:length(far_end));
w_0 = zeros(1, 128);
N = length(far_end);
M = 128;
mu = 0.25;
w = w_0;
e = zeros(1, N);
eps = 1e-6; 

for n = M:N
    x_n = fliplr(far_end(n-M+1:n));
    y_n = dot(w, x_n);
    e(n) = echo(n) - y_n;
    w = w + mu * e(n) * x_n / (dot(x_n, x_n) + eps); 
end

figure;
plot(far_end);
title('Far-End Signal');

figure;
plot(echo);
title('Echo Signal');

figure;
plot(e);
title('Error Signal');

[h_path, w_path] = freqz(path_mat.path(1,:), 1, 8000);
[h_est, w_est] = freqz(w, 1, 8000);

figure;
plot(w_path, 20*log10(abs(h_path)), 'b', 'LineWidth', 1.5);
hold on;
plot(w_est, 20*log10(abs(h_est)), 'r', 'LineWidth', 1.5);
title('Echo Path and Estimate');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('Given FIR System', 'Estimated FIR Channel');

figure;

subplot(1, 2, 1);
plot(path_mat.path(1,:), 'b', 'LineWidth', 1.5);
title('Impulse Response');
legend('Impulse Response');

subplot(1, 2, 2);
plot(w, 'r', 'LineWidth', 1.5);
title('Estimated Impulse Response');
legend('Estimated Impulse Response');
disp(['Final weights vector (w): ', num2str(w)]);