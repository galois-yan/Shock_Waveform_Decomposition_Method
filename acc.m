classdef acc < timeseries
    %
    %ACC A collection of methods dealing with acceleration time history
    %signal of mechanical shocks. This class is written as a subclass of
    %'timeseries' class, so all 'timeseries' methods can also be used.
    %
    %ACC Properties:
    %   Sf - Sample rate
    %   Time - Time column
    %   Data - Measured acceleration data column
    %   Length - Length of time series
    %
    %ACC Methods:
    %   acc - Constructor method to creat a ACC object.
    %   resample1 - Resample a time series.
    %   bandpass - A bandpass filter.
    %   plot - An overload plot function for ACC object.
    %   fft - An overload fast fourier transform function for ACC object.
    %   fit - Shock waveform decomposition method.
    %   cwt - An overload continues wavelet transform plot.
    %   dwt - An overload discrete wavelet transform plot.
    %   cumtraapz - Overload numerical integration.
    %   diff - Overload numerical difference.
    %   extend - Extending the time series for a certain period.
    %
    properties
        Sf  % Sample rate
    end
    
    methods
        function Acc = acc(Sign)
            %Acc = ACC(Sign) construct an ACC object from the measurement
            %matrix 'Sign'. The first column of 'Sign' shall be time
            %column, and the rest columns shall be acceleration data.
            
            Acc@timeseries(Sign(:,2:end), Sign(:,1)-Sign(1,1));
            Acc.Sf = (size(Sign,1)-1)/(Sign(end,1)-Sign(1,1));
            %Acc=setuniformtime(Acc,'StartTime',Sign(1,1),'EndTime',Sign(end,1));
        end
        %------------------------------------------------------------------
        function Acc = resample1(Acc, nSamples)
            %Acc = RESAMPLE1(Acc, nSamples) resample the ACC object to a
            %'nSamples' samples time series, where 'nSamples' is a scalar.
            
            t=linspace(Acc.Time(1),Acc.Time(end),nSamples)';
            y=Acc.Data;
            while 1
                try
                    y=resample(y, nSamples, length(y));
                    break;
                catch
                    y=y(1:2:end);
                end
            end
            Name=Acc.Name;
            Acc=acc([t,y]);
            Acc.Name=Name;
        end
        
        %------------------------------------------------------------------
        function Acc=bandpass(Acc, low_f, high_f)
            %Acc=BANDPASS(Acc, low_f, high_f) filters the time series
            %between [low_f, high_f] range. 'low_f' and 'high_f' are both
            %scalar.
            
            SaRa=Acc.Sf;
            
            N=500;
            %N=round(30*(SaRa/2)/(high_f-low_f));
            win=hamming(N+1);
            
            switch nargin
                case 2
                    b=fir1(N, low_f/(SaRa/2), 'low', win, 'scale');
                case 3
                    b=fir1(N,[low_f/(SaRa/2),high_f/(SaRa/2)],'bandpass', win, 'scale');
            end
            
            FilY=filter(b,1,Acc.Data);
            
            delay = round(mean(grpdelay(b)));
            t=Acc.Time(1:end-delay);
            FilY(1:delay)=[];
            Acc=acc([t, FilY]);
        end
        
        %------------------------------------------------------------------
        function plot(Acc)
            figure;
            plot(Acc.Time, Acc.Data);
            
            xlabel({'Time, s'});
            ylabel({'Acceleration, m/s^2'});
            set(gca,'GridLineStyle','--','XGrid','on','YGrid','on');
        end
        %------------------------------------------------------------------
        function [f, P1, P2, Y]=fft(Acc)
            L=length(Acc.Data);
            if rem(L,2)==0
                Y=fft(Acc.Data);
            else
                Y=fft(Acc.Data(1:end-1));
                L=L-1;
            end
            
            
            P2 = abs(Y/L);
            P1 = P2(1:L/2+1,:);
            P1(2:end-1) = 2*P1(2:end-1);
            f = Acc.Sf*(0:(L/2))/L;
            if nargout==0
                figure;
                plot(f,P1);
                title('Single-Sided Amplitude Spectrum');
                xlabel('\xi (Hz)');
                ylabel('$$|\hat{f}(\xi)|$$','Interpreter','Latex');
            end
        end
        %------------------------------------------------------------------
        function SWD = fit(Acc, varargin)
            %
            %SWD=FIT(Acc, Name, Value) retures the object of shock waveform
            %decomposition results. Optional name-value pair arguments can
            %be added.
            %
            %--------------------------------------------------------------
            %Name-Value Pairs:
            %   'FreSpace' - Spacing of frequenciey start points;
            %   2 (default) | Positive scalar
            %
            %   'IniTim' - Start point of initial time.
            %   0 (default) | Row vector
            %
            %   'TauLoc' - Start points of peak times.
            %   [0.6,1,1.4] (default) | Row vector
            %
            %   'TwoWay' - Whether consider the situation with negative
            %   peak time. '0' for 'No', and '1' for 'Yes'.
            %   0 (default) | 1
            %
            %   'ErrTol' - The tolerance of error energy ratio.
            %   0.1 (default) | Scalar between (0, 1)
            %
            %   'MinSW' - Minimum shock wavefrom components
            %   6 (default) | Positive scalar
            %
            %   'XiList' - Start points of damping ratios.
            %   logspace(-2,1,4) (default) | row vector
            %
            %   'PhiNum' - How many start points considered for phase.
            %   2 (default) | Positive scalar
            
            p=inputParser;
            addParameter(p,'FreSpace', 2);
            addParameter(p,'TauLoc', [0.6,1,1.4]);
            addParameter(p,'IniTim', 0);
            addParameter(p,'TwoWay', 0);
            addParameter(p,'ErrTol', 0.1);
            addParameter(p,'MinSW', 6);
            addParameter(p,'XiList', logspace(-2,1,4));
            addParameter(p,'PhiNum', 2);
            parse(p,varargin{:});
            
            disp(p.Results);
            
            Acc.Time=Acc.Time-Acc.Time(1);
            Acc.Data=Acc.Data-Acc.mean;
            t=Acc.Time;
            r(:,1)=Acc.Data;
            w=zeros(size(r)); % Shock waveform components.
            fval0=sum(r.^2); % Signal energy.
            
            % Set up Lower and Upper bound.
            [f, P1]=Acc.fft;
            f_U=max(f);
            tau_U=max(t); %Upper bound of time.
            [amp_U, Itau]=max(abs(Acc.Data));
            
            Lower=[-2*amp_U, 0.001,  -tau_U*p.Results.TwoWay, -tau_U*p.Results.TwoWay,  0.0001,-Inf];
            Upper=[ 2*amp_U, f_U*2*pi,2*tau_U,tau_U, 100, Inf];
            
            %%Set up multi start points.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            locs=find(P1>max(P1)/10);
            omega=2*pi*2.^(log2( f(min(locs)+1) ) : log2(p.Results.FreSpace) : ceil(log2( f(max(locs)-1) )));
            t0=p.Results.IniTim;
            tau=t(Itau)*p.Results.TauLoc;
            
            xi=p.Results.XiList;
            phi=linspace(0, 2*pi, p.Results.PhiNum+1);
            phi(end)=[];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for i=1:1000
                disp(i);
                [amp_SP,I]=max(abs(r(:,i)));
                tau_SP=t(I);
                [f,P1]=acc([t,r(:,i)]).fft;
                P1(f==0)=[];
                f(f==0)=[];
                [~,I]=max(P1);
                f_SP=f(I);
                StartPoint=[amp_SP, f_SP*2*pi, 0, tau_SP, 0.05, 0];
                
                amp=amp_SP;
                [AMP, OMEGA, T0, TAU, XI, PHI]=ndgrid(amp,omega,t0,tau,xi,phi);
                ptmatrix=[AMP(:),OMEGA(:),T0(:),TAU(:),XI(:),PHI(:)];
                
                disp([num2str(size(ptmatrix,1)) ,' local solver will run in each interation.']);
                tpoints = CustomStartPointSet(ptmatrix);
                
                problem = createOptimProblem('lsqcurvefit','x0',StartPoint,'objective',@acc.objective,...
                    'lb',Lower,'ub',Upper,'xdata',t,'ydata',r(:,i));
                ms = MultiStart('PlotFcns',@gsplotbestf);
                %ms = MultiStart('UseParallel',true,'PlotFcns',@gsplotbestf);
                
                [xmulti(i,:),fval(i)] = run(ms,problem,tpoints);
                disp((fval0-fval)./fval0);
                
                w(:,i)=acc.objective(xmulti(i,:),t);
                r(:,i+1)=r(:,i)-w(:,i);
                
                if fval(i)< p.Results.ErrTol*fval0 && i>=p.Results.MinSW
                    break
                end
                
                pause(1);
            end
            
            [~,I]=sort(sum(w.^2),'descend');
            SWD=swd(xmulti(I,:),t,r(:,1));
            
        end
        
        function cwt(Acc, cmin)
            %
            %CWT(Acc, cmin) plots the continues wavelet transform. 'cmin' is
            %to control the minimum of colormap limit, which shall be in
            %linear (normal) scale.
            
            x=Acc.Data;
            t=Acc.Time;
            [cfs,f] = cwt(x,Acc.Sf,'amor','FrequencyLimits',[100 25600]);
            %[cfs,f] = cwt(x,Acc.Sf,'amor');
            figure;
            hp = pcolor(t,f,log10(abs(cfs))); hp.EdgeColor = 'none';
            set(gca,'YScale','log');
            
            if exist('cmin','var')
                cLimit=caxis;
                cLimit(1)=log10(cmin);
                caxis(cLimit);
            end
            
            colormap jet;
            c=colorbar;
            c.TickLabels=round(10.^c.Ticks);
            
            xlabel('Time (s)'); ylabel('Frequency (Hz)');
        end
        
        function Acc=dwt( Acc0, order)
            if nargin<2
                wtecg=modwt(Acc0.Data);
                order=size(wtecg,1)-1;
            else
                wtecg=modwt(Acc0.Data,order);
            end
            mra = modwtmra(wtecg);
            mra=mra';
            
            figure;
            subplot( ceil((order+1)/2), 2, order+1);
            plot(Acc0.Time, mra(:,end));
            title('Approximation Coefficients');
            for i=1:order
                subplot(ceil((order+1)/2), 2,i)
                plot(Acc0.Time, mra(:,i));
                title(['Level ', num2str(i) ,' Detail Coefficients'])
            end
            
            Acc=Acc0;
            Acc.Data=mra;
        end
        
        function Acc=cumtrapz(Acc)
            Acc.Data=cumtrapz(Acc.Time,Acc.Data);
            if nargout==0
                Acc.plot;
            end
        end
        
        function Acc=diff(Acc)
            t=Acc.Time(1:end-1);
            y=diff(Acc.Data)*Acc.Sf;
            Acc=acc([t,y]);
            Acc.plot;
        end
        
        function Acc=extend(Acc0, duration)
            %Acc=EXTEND(Acc0, duration) extends the time history for
            %'duration' seconds.
            
            y=Acc0.Data;
            t=(0:1/Acc0.Sf:duration)';
            y=[y;zeros(length(t)-length(y),1)];
            Acc=acc([t,y]);
        end
        
    end
    
    methods (Static)
        function [ z ] = objective( x, t)
            %OBJECTIVE construct shock waveform with parameter x.
            %  'x' is the input parameter matrix, where each row stands for
            %  each waveform.
            
            Ip=x(:,4)>0;
            In=x(:,4)<0;
            Iz=x(:,4)==0;
            
            amp=x(:,1)';
            omega=x(:,2)';
            t0=x(:,3)';
            tau=x(:,4)';
            zeta=x(:,5)';
            phi=x(:,6)';
            
            t=t-t0;
            
            if sum(Ip)~=0
                z(:,Ip)=amp(Ip).*exp(zeta(Ip).*omega(Ip).*(tau(Ip)-t(:,Ip))+zeta(Ip).*tau(Ip).*omega(Ip).*(log(t(:,Ip))-log(tau(Ip)))).*cos(omega(Ip).*t(:,Ip)+phi(Ip)).*(t(:,Ip)>=0);
            end
            if sum(In)~=0
                z(:,In)=amp(In).*exp(zeta(In).*omega(In).*(t(:,In)-tau(In))+zeta(In).*tau(In).*omega(In).*(log(tau(In))-log(t(:,In)))).*cos(omega(In).*t(:,In)+phi(In)).*(1-(t(:,In)>=0));
            end
            if sum(Iz)~=0
                z(:,Iz)=amp(Iz).*exp(-zeta(Iz).*omega(Iz).*t(:,Iz)).*cos(omega(Iz).*t(:,Iz)+phi(Iz)).*(t(:,Iz)>=0);
            end
            
        end
        
    end
end
