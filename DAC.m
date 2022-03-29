% calculate wavelength (nm) from voltage (mV)
function lambda=DAC(V)
lambda=300+V/10;
end