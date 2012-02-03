function psd = psdcore(x,  ntaper, ndecimate)
%  function psd = psdcore(x,  ntaper)
%  Compute a spectral estimate of the power spectral density
%  (psd) for the time series x using sine multitapers.
%  Normalised to sampling interavl of 1.
%  
%  ntaper gives the number of tapers to be used at each frequency:
%  if ntaper is a scalar, use same value at all freqs; if a
%  vector, use ntaper(j) sine tapers at frequency j. 
%  If series length is n, psd is found at 1 + n/2 evenly spaced freqs;
%  if n is odd, x is truncted by 1.

%  ndecimate: number of psds actually computed = (1+n/2)/ndecimate;
%  these values are linearly interpolated into psd.

persistent fftz n nhalf varx

%  When ntaper is a scalar, initialize
if length(ntaper) == 1

  n = length(x(:));
  n = n - mod(n,2);    % Force series to be even in length
  nhalf = n/2;
  varx = var(x(1:n));
  ntap = ones(nhalf+1,1)*ntaper;  % Make a vector from scalar value

%  Remove mean; pad with zeros
  z = [reshape(x(1:n),[],1)  - mean(x(1:n)); zeros(n,1)];

%  Take double-length fft
  fftz = fft(z);    
else
  ntap = ntaper;
end

%  Decimation argument is optional
if nargin==2 ndecimate=1; end

%  Select frequencies for PSD evaluation
if  length(ntaper)  > 1 && ndecimate > 1
  nsum = cumsum(1 ./ntap); 
  tmp = nhalf * (nsum-nsum(1))/(nsum(end)-nsum(1));
  f =[ round(interp1(tmp,[0:nhalf]', [0 :ndecimate: nhalf]')); nhalf];
  iuniq = 1 + find(diff(f) > 0);
  f = [0; f(iuniq)];   %  Remove repeat frequencies in the list
else
  f= [0 : nhalf]';
end

%  Calculate the psd by averaging over tapered estimates
nfreq = length(f);
psd = zeros(nfreq,1);

%  Loop over frequency
for j = 1: nfreq

   m = f(j);
   tapers = ntap(m+1);
%  Sum over taper indexes; weighting tapers parabolically
   k = 1 : tapers;
   w = (tapers^2 - (k-1).^2)*(1.5/(tapers*(tapers-0.25)*(tapers+1)));
   j1=mod(2*m+2*n-k, 2*n);
   j2=mod(2*m+k, 2*n);
   psd(j) = w * abs(fftz(j1+1)-fftz(j2+1)).^2;

end

%  Interpolate if necessary to uniform freq sampling
if length(ntaper) > 1 && ndecimate > 1
  psd = interp1(f,psd, [0:nhalf]');
  end

%  Normalize by variance
area = (sum(psd) - psd(1)/2 - psd(end)/2)/nhalf;  % 2*Trapezoid
psd = (2*varx/area)*psd;


