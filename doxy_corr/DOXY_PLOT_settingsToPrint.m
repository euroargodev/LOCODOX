% DOXY_PLOT_settingsToPrint prepares and print figures.
%
% SYNTAX
% [hFig] = DOXY_PLOT_settingsToPrint(hFig,Work,saveFile)
%
% DESCRIPTION
% DOXY_PLOT_settingsToPrint prepares the figure to be printed, by settings
% the dimension, making it visible "off", then prints it and releases in the
% initial state before displays the figure again.
%
% INPUT
%   hFig                handles of the figure to be printed
%
%   Work (struct)      float working structure, issued and computed from
%                      the argo data and the reference data (WOA or REF).
%                      Carries plot informations among other things.
%                      Example:
%                      Work = 
% 
%                                  readme: [1x1 struct]
%                                    unit: 'mumol/kg'
%                              R2threshold: 0.8000
%                                     wmo: 1901205
%                                  sensor: 'Optode'
%                                makePlot: 1
%                                savePlot: 1
%                                 dirPlot: '/home1/homedir4/perso/ebrion/NAOS/NAOS_2015(2)/plots/REF'
%                           pres_adjusted: [1x1 struct]
%                           temp_adjusted: [1x1 struct]
%                           psal_adjusted: [1x1 struct]
%                                DOXY_WOA: [1x1 struct]
%
%   saveFile           full name of the figure when printed.
%
% OUTPUT
%   hFig                handles of the figure to be printed
%
% CALL :
%
% SEE ALSO
%   

% HISTORY
%   $created: 25/07/2018 $author: Emilie Brion, Altran Ouest
%   $Revision: version $Date: $author:

function [hFig] = DOXY_PLOT_settingsToPrint(hFig,Work,saveFile)


% =========================================================================  
%% Initialisation : get the initial parameters of the figure
% =========================================================================  

%set(hFig, 'PaperPositionMode','manual')
set(hFig,'Units','centimeters');
set(hFig,'PaperUnits','centimeters');
initpos = get(hFig,'Position');
initPaperPos = get(hFig,'PaperPosition');
initPaperSize = get(hFig,'PaperSize');

% =========================================================================  
%% Modify the figure to make it printable
% =========================================================================  
set(hFig,'Units','normalized');
set(hFig,'Visible','off');
set(hFig,'Position',[0 0 0.99 0.99]);
set(hFig,'Units','centimeters');
set(hFig,'PaperUnits','centimeters');
finalpos = get(hFig,'Position');
set(hFig,'PaperPosition',[0 0 finalpos(3:4)]);
set(hFig,'PaperSize',finalpos(3:4))
set(hFig, 'PaperPositionMode', 'auto') %02/07/19 marine

% =========================================================================  
%% Print the figure
% =========================================================================
if length(Work.formattype)>1
    for i = 1:length(Work.formattype)
        if isempty(strfind(saveFile,Work.formattype{i}))
            saveFile=strcat(saveFile,strrep(Work.formattype{i},'-d','.'));
        end
        print(hFig,Work.formattype{i},Work.resol,saveFile);
        if Work.savePlotFig
            suf=strrep(Work.formattype{i},'-d','.');
            saveFileFig=strrep(saveFile,suf,'.fig');
            % DOXY_PLOT_data_corr_*.fig 3x3
            % DOXY_PLOT_interpolation*.fig ==> time drift 2 x 1
            saveifig=0;
            if ~isempty(strfind(saveFileFig,'DOXY_PLOT_corr'))
                saveifig=1;
            end
            if ~isempty(strfind(saveFileFig,'DOXY_PLOT_data_corr'))
                saveifig=1;
            end
            if ~isempty(strfind(saveFileFig,'DOXY_drift_on'))
                saveifig=1;
            end
            if ~isempty(strfind(saveFileFig,'DOXY_PLOT_interpolation'))
                saveifig=1;
            end
            if saveifig
                savefig(saveFileFig);
            end
        end
    end
else
    if isempty(strfind(saveFile,Work.formattype{1}))
        saveFile=strcat(saveFile,strrep(Work.formattype{1},'-d','.'));
    end
    print(hFig,Work.formattype{1},Work.resol,saveFile);
    if Work.savePlotFig
        suf=strrep(Work.formattype{1},'-d','.');
        saveFileFig=strrep(saveFile,suf,'.fig');
        saveifig=0;
        if ~isempty(strfind(saveFileFig,'DOXY_PLOT_corr'))
            saveifig=1;
        end
        if ~isempty(strfind(saveFileFig,'DOXY_PLOT_data_corr'))
            saveifig=1;
        end
        if ~isempty(strfind(saveFileFig,'DOXY_drift_on'))
            saveifig=1;
        end
        if ~isempty(strfind(saveFileFig,'DOXY_PLOT_interpolation'))
            saveifig=1;
        end
        if saveifig
            savefig(saveFileFig);
        end
    end
end

% =========================================================================  
%% setting back the figure
% =========================================================================  
set(hFig,'PaperPosition',initPaperPos);
set(hFig,'PaperSize', initPaperSize);
set(hFig,'Position',initpos);
set(hFig,'Visible','on');