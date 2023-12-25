function myPrint(FileName,width,height,FontSize,LineWidth,MarkerSize,LegendFontSize,Res)

    % Set default values
    if ~exist('MarkerSize','var') || isempty(MarkerSize)
        MarkerSize = 2;
    end
    if ~exist('LineWidth','var') || isempty(LineWidth)
        LineWidth = 0.9;
    end
    if ~exist('FontSize','var') || isempty(FontSize)
        FontSize = 8;
    end
    if ~exist('LegendFontSize','var') || isempty(LegendFontSize)
        LegendFontSize = 6;
    end
    if ~exist('Res','var') || isempty(Res)
        Res = '-r1200';
    end

    % Get all currently open figures
    figHandles = findall(groot, 'Type', 'figure');

    for h = figHandles'
        figure(h)
        Axis_list = findall(h.Children,'Type','Axes');
        for jj = 1:length(Axis_list)
            Axis_list(jj).FontSize = FontSize;
            Axis_list(jj).LineWidth = 0.25;
            try
                [Axis_list(jj).Children.LineWidth] = deal(LineWidth);
                [Axis_list(jj).Children.MarkerSize] = deal(MarkerSize);
            catch
            end
            Axis_list(jj).FontName = 'TimesNewRoman';
            for ii = 1:length(Axis_list(jj).YAxis)
                [Axis_list(jj).YAxis(ii).Label.Interpreter] = deal('latex');
                [Axis_list(jj).YAxis(ii).TickLabelInterpreter] = deal('latex');
            end
            [Axis_list(jj).XAxis.TickLabelInterpreter] = deal('latex');
            [Axis_list(jj).XAxis.Label.Interpreter] = deal('latex');
            try
                [Axis_list(jj).ZAxis.TickLabelInterpreter] = deal('latex');
                [Axis_list(jj).ZAxis.Label.Interpreter] = deal('latex');
            catch
            end
            
            Axis_list(jj).Title.Interpreter = 'latex';
%             try
%                 LegendStrings = Axis_list(jj).Legend.String;
%                 delete(Axis_list(jj).Legend);
%                 axes(Axis_list(jj));
%                 legend(LegendStrings,'Interpreter','latex','FontSize',LegendFontSize,'Location','best');
%             catch
%             end
            
            ColorBar = findall(h.Children,'Type','ColorBar');
            if ~isempty(ColorBar)
                ColorBar.TickLabelInterpreter = 'latex';
            end
        end
        
        h.Clipping = false;
        h.PaperUnits = 'inches';
        h.PaperPosition = [0,0,width,height];
        h.PaperSize = [width,height];
        h.Units = 'inches';
        h.Position = [0,0,width,height];
%         h.InnerPosition =[0,0,width,height];
        % h.Renderer = 'painters';

%         print([FileName,'_',num2str(h.Number)],'-dpdf');
%         print(Res,[FileName,'_',num2str(h.Number)],'-dpng');
%         savefig(h,[FileName,'_',num2str(h.Number)]);
        
    end

end

