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
%           30.01.2025     fingure print + debug filename by Virginie Racape

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
% set(hFig,'Units','normalized'); % Jan26, vr - pokapok, commented
% set(hFig,'Visible','off');
% set(hFig,'Position',[0 0 0.99 0.99]);
% set(hFig,'Units','centimeters');
% set(hFig,'PaperUnits','centimeters');
% finalpos = get(hFig,'Position');
% set(hFig,'PaperPosition',[0 0 finalpos(3:4)]);
% set(hFig,'PaperSize',finalpos(3:4))
% set(hFig, 'PaperPositionMode', 'auto') %02/07/19 marine 

% Jan26, vr - pokapok : change to be screen independant
set(hFig,'Visible','off');
pos = get(hFig,'Position');
figsize = pos(3:4);

%paper orientation
Lpos = find(figsize == max(figsize));
coef = 32 / figsize(Lpos(1));
finalpos = [0 0 (figsize .* coef)];
set(hFig,'Position',finalpos)
set(hFig, 'PaperType','A4');

% =========================================================================  
%% Print the figure
% =========================================================================
if length(Work.formattype)>1
    for i = 1:length(Work.formattype)
        if isempty(strfind(saveFile,Work.formattype{i}))
            saveFile_f=strcat(saveFile,strrep(Work.formattype{i},'-d','.')); % Jan26, vr - pokapok : change file parameter
        end
        %figure(hFig);
        print(hFig,Work.formattype{i},Work.resol,saveFile_f); % Jan26, vr - pokapok : change file parameter
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
                figure(hFig);
                savefig(saveFileFig);
            end
        end
    end
else
    if isempty(strfind(saveFile,Work.formattype{1}))
        saveFile_f=strcat(saveFile,strrep(Work.formattype{1},'-d','.')); % Jan26, vr - pokapok : change file parameter
    end
    %figure(hFig);
    print(hFig,Work.formattype{1},Work.resol,saveFile_f); % Jan26, vr - pokapok : change file parameter
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
            figure(hFig);
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