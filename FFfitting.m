% This function takes mean radius, standard deviation in radii, background
% and a scaling factor as param vector, and the Q of inteerst and returns a
% theoretical form factor. The distribution of particle sizes is assumed to
% be gaussian.
function yfit = FFfitting(param,Q)
rMean = param(1);
rSpread = param(2);
ptNumber = 200;

d_r = rSpread*8/ptNumber;
rRange = rMean-4*rSpread:d_r:rMean+4*rSpread;
f_r = normpdf(rRange,rMean,rSpread); % Generate normal distribution


bkg = param(3);
scale = param(4);
IQ=zeros(1,length(Q));

for k = 1:length(Q);
    
    QR = Q(k)*rRange;
    FQ = 3*(sin(QR)-QR.*cos(QR))./QR.^3;
    IQ(k) = sum(f_r.*(rRange.^6).*FQ.^2*d_r);
    
end

yfit = scale*IQ.' + bkg;