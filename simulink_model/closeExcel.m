function closeExcel()
    try
        app = actxGetRunningServer('Excel.Application');
        app.DisplayAlerts=false;
        if app.Workbooks.Count>0
            app.ActiveWorkbook.Save();
            fprintf('Final iteration results saved to Excel.\n');
        end
        app.Quit();
        delete(app);
        fprintf('Excel shut down successfully.\n');
    catch
    end
end