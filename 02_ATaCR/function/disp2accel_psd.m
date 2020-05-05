function psd = disp2accel_psd(DISP,NFFT,dt)
samprate = 1/dt;
VEL= time_differentiate(DISP, samprate);
ACC= time_differentiate(VEL, samprate);
xdft = fft(ACC,NFFT);
xdft = xdft(1:NFFT/2+1);
psdx = (1/(samprate*NFFT)).*abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
psd = smooth(psdx,100);

end

function out = time_differentiate(in,sps)
out = zeros(1,length(in));
for i=1:length(in) - 1
    out(i) = (in(i+1) - in(i)) * sps;
end
out(length(in)) = out(length(in) - 1);

end