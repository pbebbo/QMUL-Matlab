function plotRegressionResults(x, y, data_pts, d, y_pred, x_RMS, y_RMS)
    % Plot regression results in subplots
    subplot(2,2,1);
    plot(x, y_pred.under, 'b-', x, y, 'g-', data_pts(1:d,1), data_pts(1:d,2), 'ro');
    axis([0,1,-1.5,1.5]);
    title('d=2 Underfitting');
    
    subplot(2,2,2);
    plot(x, y_pred.close, 'b-', x, y, 'g-', data_pts(1:d,1), data_pts(1:d,2), 'ro');
    axis([0,1,-1.5,1.5]);
    title('d=3 close fitting');
    
    subplot(2,2,3);
    plot(x, y_pred.over, 'b-', x, y, 'g-', data_pts(1:d,1), data_pts(1:d,2), 'ro');
    axis([0,1,-1.5,1.5]);
    title('d=9 overfitting');
    
    subplot(2,2,4);
    plotRMSErrors(x_RMS, y_RMS);
end 