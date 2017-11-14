% This file defines the data class for SAXS data.
classdef SAXSdata
    properties
        SN;
        Q;
        Intensity;
        Error;
        
        bSN;
        bQ;
        bIntensity;
        bError;
        sf;
    end
    methods
        function sxd = SAXSdata(t,b,scaling_factor)
            if nargin >0
                s = strcat('*',num2str(t,'%05i'),'_00001.dat');
                files = dir(s);
                if isempty(files)== 0
                    ALL_DATA=importdata(files.name,'\t',13);
                    data=ALL_DATA.data;
                end
                if isempty(files) == 1
                    s = strcat('*',num2str(t,'%05i'),'.avg');
                    files = dir(s);
                    ALL_DATA=importdata(files.name,'\t',8);
                    data=ALL_DATA.data;
                end
                
                if isempty(files) == 1
                    error('Data file with .dat extension does not exist in this folder')
                end
                
                
                sxd.Q = data(:,1);
                sxd.Intensity = data(:,2);
                sxd.Error = data(:,3);
                sxd.SN=t;
                sxd.sf=scaling_factor;
                if b~= 0
                    s = strcat('*',num2str(b,'%05i'),'_00001.dat');
                    files = dir(s);
                    if isempty(files)== 0
                        ALL_DATA=importdata(files.name,'\t',13);
                        data=ALL_DATA.data;
                    end
                    if isempty(files) == 1
                        s = strcat('*',num2str(b,'%05i'),'.avg');
                        files = dir(s);
                        ALL_DATA=importdata(files.name,'\t',8);
                        data=ALL_DATA.data;
                    end
                    
                    if isempty(files) == 1
                        error('Data file with .dat extension does not exist in this folder')
                    end
                    
                    sxd.bQ = data(:,1);
                    sxd.bIntensity = data(:,2);
                    sxd.bError = data(:,3);
                    sxd.bSN=b;
                else
                    sxd.bSN=0;
                end
                
            end
            
            
        end
        function background_subtracted = subtract(sample)
            if sample.bSN==0
                background_subtracted=sample;
            else
                background_subtracted = SAXSdata;
                background_subtracted.Q = sample.Q;
                background_subtracted.Intensity = sample.Intensity - sample.sf * sample.bIntensity;
                background_subtracted.Error = sqrt(sample.Error.^2 + (sample.sf.^2*sample.bError.^2));
            end
        end
        function fi = plotraw(sample)
            %j =  loglog(sample.Q, sample.Intensity, '-*');
            j =  ploterr(sample.Q, sample.Intensity, [],sample.Error,'-*','logxy');
            fi = ancestor(j(1), 'figure');
            set(j,'LineWidth',2)
            axis([0.01 0.5 0 inf])
            axis square
            xlabel('q ($\rm{\AA} ^{-1}$ )','interpreter','LaTex','FontSize',30)
            ylabel('Intensity (arb. units)','FontSize',30)
            ax=gca;
            set(ax,'XScale','log')
            set(ax,'YScale','log')
            ax.XAxis.FontName = 'Times New Roman';
            ax.YAxis.FontName = 'Times New Roman';
            ax.XAxis.FontSize=20;
            ax.YAxis.FontSize=20;
            ax.Title.FontSize=20;
            ax.LineWidth =2;
            set(ax, 'box','on')
            set(findall(gcf,'type','text'),'FontName','Times New Roman','FontSize',30)
        end
        function fi = plotbs(sample)
            if sample.bSN==0
                error('Sample has no background data provided')
            else
                bs = subtract(sample);
                h = ploterr(bs.Q, bs.Intensity, [],bs.Error,'-*','logxy');
                set(h,'LineWidth',2)
                fi = ancestor(h(1), 'figure');
                axis([0.01 0.5 0 inf])
                axis square
                xlabel('q ($\rm{\AA} ^{-1}$ )','interpreter','LaTex','FontSize',30)
                ylabel('Intensity (arb. units)','FontSize',30)
                ax=gca;
                ax.XAxis.FontName = 'Times New Roman';
                ax.YAxis.FontName = 'Times New Roman';
                ax.XAxis.FontSize=20;
                ax.YAxis.FontSize=20;
                ax.Title.FontSize=20;
                ax.LineWidth =2;
                set(ax,'XScale','log')
                set(ax,'YScale','log')
                set(ax, 'box','on')
                set(findall(gcf,'type','text'),'FontName','Times New Roman','FontSize',30)
            end
        end
    end
end
