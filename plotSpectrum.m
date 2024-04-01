function []=plotSpectrum(x,fs,pc,offset)

% if isempty(pc)
%     pc='b';
% end

x0_f=fftshift(fft(x));
ff=fs/2*linspace(-1,1,length(x0_f));
x0_f_db=20*log10(abs(x0_f));
%x0_f_db=x0_f_db-max(x0_f_db);

if isempty(pc)
    plot(ff/1e6,x0_f_db+offset)
else
    plot(ff/1e6,x0_f_db+offset,'color',pc)
end
grid on;xlabel('Freq (MHz)')