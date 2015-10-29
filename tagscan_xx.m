%tagscan_xx
%Produce velocity distribution as a function of nu for different xx values.
function Data = tagscan_xx( xx,nu,cf,sepz,a,tp,B )
%**************************************************************************
wcf = 2*pi*1000*cf;                  % Angular chop frequency
wc = 9580*B/138;
N = 30000;                           % Discretize time into 30000 intervals
t = ((1:N)* 1e-2)/wc;
t2 = [-t(end:-1:2) t];
%**************************************************************************
Data = zeros(length(xx),length(nu)); % Create a 11x101 array to store data
for i = 1:length(xx)
    xi = xx(i);
    parfor j = 1:length(nu)
        nj = nu(j);
        S = ptag(xi,nj,sepz,a,tp,wc);
        Data(i,j) = abs(sum(exp(sqrt(-1)*wcf*t2).*xcorr(S,square(wcf*t),'unbiased')));
    end
    message = sprintf('%d of 11 is finished.',i);
    disp(message);
end
end

