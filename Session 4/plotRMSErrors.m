function plotRMSErrors(x_RMS, y_RMS)
    % Plot RMS errors
    plot(x_RMS, y_RMS.data_under, 'bx-', x_RMS, y_RMS.test_under, 'bx--',...
         x_RMS, y_RMS.data_close, 'rx-', x_RMS, y_RMS.test_close, 'rx--',...
         x_RMS, y_RMS.data_over, 'gx-', x_RMS, y_RMS.test_over, 'gx--');
    axis([0, max(x_RMS), 0, 1]);
    title('RMS error');
    legend('d=2 - data', 'd=2 - test', 'd=3 - data', 'd=3 - test',...
           'd=9 - data', 'd=9 - test');
end 