classdef swd
    %
    %SWD A collection of methods dealing with shock waveform components
    %decomposed from time histories.
    %
    %SWD Properties:
    %   x - Shock wavefrom parameters. Each row corresponds to a shock
    %   waveform.
    %   t - Time column of original time series
    %   y - Data column of original time series
    %   energy - The signal energy of each shock waveform component
    %   energyRatio - The energy ratio (compared with original signal) of
    %   each shock waveform component
    %   kappa - Normalized peak time
    %   n - Number of shock waveform components
    %   eta - The number of necessary shock wavefrom components to reach a
    %   90% reconstruction.
    %   w - The shock waveform matrix, where each column is the time
    %   history of a single shock waveform component.
    %
    %SWD Methods:
    %   swd - Constructor method to creat a SWD object.
    %   table - Show the shock waveform parameters in the form of a table.
    %   yhat - Reconstruct a time history from shock wavefrom components.
    %   subplot - Show all shock waveforms in the form of subplot.
    
    properties
        x
        t
        y
    end
    properties (Dependent)
        energy
        energyRatio
        kappa
        n
        eta
        w
    end
    
    methods
        function obj=swd(x,t,y)
            %
            %obj = SWD(x, t, y) constructes the SWD object from SW parameters,
            %time and data columns. The 'x' is a row-wise parameter matrix;
            %the 't' and 'y' are time and data columns, respectively.
            
            obj.x=x;
            obj.t=t;
            obj.y=y;
        end
        
        function value=get.w(obj)
            value=acc.objective(obj.x, obj.t);
        end
        
        function value=get.energy(obj)
            value=sum(obj.w.^2)';
        end
        
        function value=get.energyRatio(obj)
            value=obj.energy./sum(obj.y.^2);
        end
        
        function value=get.kappa(obj)
            value=obj.x(:,4).*(obj.x(:,2)/(2*pi));
        end
        
        function value=get.eta(obj)
            [value,~]=find(cumsum(obj.energyRatio)>=0.9, 1 );
            if isempty(value)
                value=size(obj.x, 1);
            end
        end
        
        function value=get.n(obj)
            value=obj.x(:,2).*obj.x(:,4).*obj.x(:,5)+1;
        end
        
        function T=table(obj)
            xmulti=[obj.x, obj.energyRatio, obj.kappa];
            
            varphi=rem(xmulti(:,6),2*pi);
            varphi(varphi<0)=varphi(varphi<0)+2*pi; %#ok<*PROP>
            xmulti(:,6)=varphi;
            xmulti(:,2)=xmulti(:,2)./(2*pi);
            xmulti(:,3)=xmulti(:,3)*1000;
            xmulti(:,4)=xmulti(:,4)*1000;
            xmulti(:,7)=xmulti(:,7)*100;
            
            xmulti=round(xmulti,3);
            T=array2table(xmulti,'VariableNames',{'Amplitude','Frequency_Hz',...
                't0_ms','tau_ms','zeta','varphi','epsilon','kappa',});
            
        end
        
        function Acc=yhat(obj, I)
            %Acc=YHAT(obj, I) reconstruct a time history from shock
            %wavefrom components. 'I' is the index vector of selected
            %waveforms. The reconstructed time series is shown as a 'acc'
            %object.
            
            Acc=acc([obj.t,sum(obj.w,2)]);
            
            if nargin == 2
                %wSelect=obj.w(:,I);
                Acc=acc([obj.t,sum(obj.w(:,I),2)]);
            end
            
        end
        
        function subplot(obj, pages)
            %SUBPLOT(obj, pages) show shock waveforms in the form of
            %subplot. 'pages' define how many figure is going to show. Each
            %figure can show 10 shock wavefrom components.
            
            if exist('pages', 'var')
                pages=ceil( size(obj.w,2)/10);
            else
                pages=1;
            end
            
            for j=fliplr(1:pages)
                figure('Units','normalized','OuterPosition',[0.25,0.05,0.4,0.95]);
                epsilon=obj.energyRatio;
                for i=( 1+(j-1)*10 : min(10*j, size(obj.w,2)) )
                    subplot(5,2, rem(i-1,10)+1);
                    plot(obj.t,obj.w(:,i));
                    xlabel({'Time, s'});
                    ylabel({'Acceleration, m/s^2'});
                    ylim([-1.1*max(abs(obj.w(:,1))),1.1*max(abs(obj.w(:,1)))]);
                    xlim([obj.t(1),obj.t(end)]);
                    legend(['\epsilon_{w',num2str(i),'}=',sprintf('%0.1f',epsilon(i)*100),'%']);
                    set(gca,'GridLineStyle','--','XGrid','on','YGrid','on');
                end
            end
            
        end
        
    end
    
end

