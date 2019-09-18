function r_y_double = ACF_two_sided(y,N)

r_y = zeros(1,N);

for k = 1:N
    if k < N 
        temp = 0;
        for n = 1:N-k-1
            temp = temp + 1/N * y(n + k - 1)*y(n); 
        end
        r_y(k) = temp;
    else
       r_y(k) = 0;
    end
end

% plot(t_axis,r_y); title('Estimation of ACF, 1st order Butterworth'); xlabel('k') ylabel('Correlation')

% för att få dubbelsidig ACF
r_y_half = r_y(2:N);                        % length N-1
r_y_half = flip(r_y_half);
r_y_double = zeros(1,2*N-1);                % length 2*N-1
r_y_double(1:N-1) = r_y_half; 
r_y_double(N:2*N-1) = r_y;