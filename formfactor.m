classdef formfactor < SAXSdata
    properties
        fQ;
        fIntensity;
        RD;
        fval;
        residual;
        exitflag;
        output;
    end
    methods
        
        function yfit = formfactor(varargin)
            yfit=yfit@SAXSdata(varargin{:});
            data=subtract(yfit);
            outliers = excludedata(data.Q,data.Intensity,'domain',[0.0075 0.5]);
            x1=data.Q(~outliers);
            y1=data.Intensity(~outliers);
            yErr1=data.Error(~outliers);
            % fitting to size Gaussian Distribution
            
            %         rMean    rSpread    backgroud       scale
            
            start=   [45;      5;         1e-4;           4e-11];
            LB=      [25;      2;         1e-8;           1e-15];
            UB=      [75;      6;         1e-2;           1e-6];
            
            %% Fitting the data to scattering from a gaussian size distribution (form factor only)
            options=optimset('MaxIter',5000,'TolFun',1e-16,'Display','on','TolX',1e-16,'MaxFunEvals',20000);
            [RD,fval,residual,exitflag,output]=lsqcurvefit(@(param,Q) FFfitting(param,Q)./yErr1,start, x1,(y1./yErr1),LB,UB,options);
            yfit.fIntensity = FFfitting(RD,x1);
            yfit.fQ = x1;
            yfit.RD = RD;
            yfit.fval =fval;
            yfit.residual =residual;
            yfit.exitflag = exitflag;
            yfit.output=output;
        end
        
        function fi = plotFF(sample)
            j =  plot(sample.Q, sample.Intensity,'-*',sample.fQ, sample.fIntensity,'LineWidth',2);
            fi = ancestor(j, 'figure');
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
            set(ax, 'box','on')
            axis square
            set(ax,'YScale','log')
            set(ax,'XScale','log')
            legend('Raw data','Fit form factor','FontWeight','normal')
            set(legend, 'FontSize',18)
            set(findall(gcf,'type','text'),'FontName','Times New Roman')
        end
    end
end