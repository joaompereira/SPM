function [] = pdfprint(filename,figure)
% Print a figure in pdf (for nice vectorized pictures)
% Solution in https://www.mathworks.com/matlabcentral/answers/12987-how-to-save-a-matlab-graphic-in-a-right-size-pdf

    if ~exist('figure', 'var') || isempty(figure)
        figure = gcf;
    end
    
    set(figure ,'Units','Inches');
    pos = get(figure ,'Position');
    set(figure ,'PaperPositionMode','Auto',...
                'PaperUnits','Inches',...
                'PaperSize',[pos(3), pos(4)]);
    print(figure, filename,'-dpdf','-r0');




end