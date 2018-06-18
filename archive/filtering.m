function [H]=filtering(filter,nfp);

 d=0.95;
 filt=[2*(0:nfp/2)./nfp]';
 w = [2*pi*(0:length(filt)-1)/nfp]';   % frequency axis up to Nyquist
% 
%filter
    switch filter
        case 'none'
            filt=[ones(1,nfp/2+1)]'; % thst is : Do nothing
             'projection filter: none '
        case 'ram-lak'
            filt=[2*(0:nfp/2)./nfp]';
            'projection filter: ram-lak '
        case 'shepp-logan'
            % be careful not to divide by 0:
            filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2*d))./(w(2:end)/(2*d)));
             'projection filter: shepp-logan '
        case 'cosine'
           
            filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*d));
        case 'hamming'
            filt(2:end) = filt(2:end) .* (.54 + .46 * cos(w(2:end)/d));
        case 'hann'
            filt(2:end) = filt(2:end) .*(1+cos(w(2:end)./d)) / 2;
             'projection filter: hann '
        case 'blackman'
            filt(2:end) = filt(2:end) .*(0.42 - 0.5*cos(w(2:end)./d) + ...
                0.08*cos(2.*w(2:end)./d));
             'projection filter: blackman '
        case 'welch'
            filt(2:end) = filt(2:end) .*(1-(w(2:end)/d).^2);
            'projection filter : welch'
        case 'connes'    
            filt(2:end) = filt(2:end) .*(1-(w(2:end)/d).^2).^2;
            
            'projection filter : connes '
        case 'blackmanharris' % Blackman-Harris window
             % check !!
            t=w(2:end)./d;
            filt(2:end)= filt(2:end).*...
                (  0.35875-.48829*cos(t)+.14128*cos(2*t)-.01168*cos(3*t));  % half window
        case 'udf' % user defined 
            %filt=[ones(1,nfp/2+1)]'; % that is : Do nothing
              filt=(2*(0:nfp/2)./nfp)'; % r-l
              filt( filt < 0.2 ) = 0.2;
             'projection filter : user defined filter'
        case 'expwin'
        beta=1; N=nfp; M = (N-1)/2; 
        for k = 0:M
        n = k-M;
        filt(k+1) = exp(beta*sqrt(1-4*n*n/(N-1)/(N-1)))/exp(beta);
        end
            
        case 'parzen'
            filt=parzenwin_1(nfp);
            'projection filter: parzen '
        case 'kaiser'
            %  Joe Henning - Dec 2013
            beta=1 ;% 1 kaiser,0 Rectangular,5 Hann,6 Hamming,12.2 Blackman-Harris
            p=1;  N=nfp; M = (N-1)/2;
            w = [];
            for k = 0:M
                n = k-M;
                w(k+1) = besseli(0,beta*sqrt(1-4*n*n/(N-1)/(N-1)))/besseli(0,beta);
                % w(N-k) = w(k+1);
            end
            w = w.^p;  w = w(:)/max(w);
            'projection filter: kaiser '
        case 'Lanczos'
            n=length(filter) ;
            if ~rem(n,2)
                % Even length window
                m = n/2;
                x = (0:m-1)'/(n-1);
                w = sinc(2*x - 1);
                %  w = [w; w(end:-1:1)];
            else
                % Odd length window
                m = (n+1)/2;
                x = (0:m-1)'/(n-1);
                w = sinc(2*x - 1);
                %  w = [w; w(end-1:-1:1)];
            end
            filt=w;
            
            
        otherwise
            error('Invalid selection: the filter selected in not listed.');
    end
    
    %%filt
    
    H=[filt;filt(end-1:-1:1)];
    H(nfp+1:end)=[];
    % show frquency filter 
    % figure, plot(H,'r+'), title(['frequency filter: ', filter])

end