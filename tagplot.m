%**************************************************************************
B = 3500;          % Magnetic field strength (Gauss)
cf =  1;           % Chop frequency (kHz)
xx = -1:.2:1;      % Tag beam spatial location
nu = -1:.02:1;     % Laser detuning frequency (GHz)
sepz = 7.8;        % Separation in z direction (cm)
a = 0.1;           % Laser beam width (cm)
tp = 0.2;          % Temperature in perpendicular direction
%**************************************************************************
Data = tagscan_xx( xx,nu,cf,sepz,a,tp,B );

color = ['k--', 'm--', 'g--', 'r--', 'b--', 'k', 'b', 'r', 'g', 'm', 'k:'];
figure;
for i = 1:length(xx)
    hold on;
    plot(nu,Data(i,:),color(i));
end
xlabel('nu (GHz)');
ylabel('Amplitude');
title(['B = ' num2str(B) ' Gauss, cf = ' num2str(cf) ' kHz, z = ' num2str(sepz) ' cm, bandwidth = ' num2str(a) 'vary xx']);