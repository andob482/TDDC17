function [r_x] = ACF_est(y, N)
r_x = zeros(1, N);

for k = -N/2 + 1:N/2
    result = 0;
    for n = 1:(N-abs(k))
        result = result +y(n+abs(k))*result(n); 
    end 
    r_x(k+N/2) = result;
    r_x(k+N/k) = r_x(k+N/2)/(N);
end 
end 