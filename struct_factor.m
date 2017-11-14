classdef struct_factor < SAXSdata
    properties
        SF; % The calculated structure factor
        SFn;
        Error_SF;
        Q_SF; % Q range used for calculating the structure factor
        ff; % Calculated form factor
        Error_SFn; % Normalized error
    end
    methods
        function strf = struct_factor(param,Q,varargin)
            strf = strf@SAXSdata(varargin{:});
            strf.ff = FFfitting(param,Q);
            strf_bkg = subtract(strf);
            %% Subsetting the data
            Q_cut=strf.Q(strf.Q>=Q(1)& strf.Q<= Q(end));
            I_cut=strf_bkg.Intensity(strf.Q>=Q(1)& strf.Q<= Q(end));
            Err_cut=strf_bkg.Error(strf.Q>=Q(1)& strf.Q<= Q(end));
            strf.Q_SF = Q_cut;
            %% Raw structure factor
            yff = interp1(Q, strf. ff, Q_cut);
            strf.SF = I_cut./yff;
            strf.Error_SF=Err_cut./yff;
            
            %% Normalizing the structure factor
            %x0=0.0485; % for all except 192
            x0=0.055;    % for 192
            xf=0.4;
            nuf=sum(Q_cut < xf);
            Sf=mean(strf.SF(nuf:end));
            nu=sum(Q_cut<x0);
            slope=strf.SF(nu)/Q_cut(nu);
            x_low=Q_cut(1:nu);
            S_low=slope*x_low';
            strf.SFn=strf.SF;
            strf.SFn(1:nu)=S_low;
            strf.SFn=strf.SFn./Sf;
            strf.Error_SFn = strf.Error_SF./Sf;
        end
        function fi = plotSFraw(sample)
            j =  ploterr(sample.Q_SF, sample.SF,[],sample.Error_SF,'-*');
            fi = ancestor(j, 'figure');
%             set(j,'LineWidth',2)
            axis([0.01 0.5 0 inf])
            axis square
            xlabel('q ($\rm{\AA} ^{-1}$ )','interpreter','LaTex','FontSize',30)
            ylabel('S(q) (arb. units)','FontSize',30)
            ax=gca;
            ax.XAxis.FontName = 'Times New Roman';
            ax.YAxis.FontName = 'Times New Roman';
            ax.XAxis.FontSize=20;
            ax.YAxis.FontSize=20;
            ax.Title.FontSize=20;
            ax.LineWidth =2;
            set(ax, 'box','on')
            axis square
            set(findall(gcf,'type','text'),'FontName','Times New Roman')
        end
        function fi = plotSFn(sample)
            j =  ploterr(sample.Q_SF, sample.SFn,[],sample.Error_SFn,'*');
            fi = ancestor(j, 'figure');
%             set(j,'LineWidth',2)
            axis([0.01 0.5 0 inf])
            axis square
            xlabel('q ($\rm{\AA} ^{-1}$ )','interpreter','LaTex','FontSize',30)
            ylabel('S(q) (arb. units)','FontSize',30)
            ax=gca;
            ax.XAxis.FontName = 'Times New Roman';
            ax.YAxis.FontName = 'Times New Roman';
            ax.XAxis.FontSize=20;
            ax.YAxis.FontSize=20;
            ax.Title.FontSize=20;
            ax.LineWidth =2;
            set(ax, 'box','on')
            set(findall(gcf,'type','text'),'FontName','Times New Roman')
        end
        function gr = radial_dist(sample)
            r=0:1:200;
            dQ=sample.Q_SF(end)-sample.Q_SF(1)/length(sample.Q_SF);
            gr=zeros(1,length(r));
            np=sqrt(2)/(9^3);
            A=1/(2*(pi^2)*np);
            for i=1:length(r)
                I= sum((sample.SFn-1).* sample.Q_SF.^2.* sin(sample.Q_SF*r(i))./(sample.Q_SF*r(i))*dQ);
                gr(i)= 1+ A* I;
            end
        end
    end
end
