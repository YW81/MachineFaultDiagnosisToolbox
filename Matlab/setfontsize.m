function setfontsize(ftsz)

set(get(gca,'XLabel'),'FontSize',ftsz,'fontname', 'Times New Roman');
set(get(gca,'YLabel'),'FontSize',ftsz,'fontname', 'Times New Roman');
set(get(gca,'title'),'FontSize',ftsz,'fontname', 'Times New Roman');
set(gca,'FontSize',ftsz,'fontname', 'Times New Roman')
end
