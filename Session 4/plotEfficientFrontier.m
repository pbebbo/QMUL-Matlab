function plotEfficientFrontier(sig, r)
    % Plot the efficient frontier
    % Inputs:
    %   sig: risk values
    %   r: return values
    
    figure;
    plot(sig, r);
    xlabel('Risk, \sigma^2', 'fontsize', 18);
    ylabel('Return, R', 'fontsize', 18);
    title('Efficient Frontier', 'fontsize', 18);
    grid on;
end 