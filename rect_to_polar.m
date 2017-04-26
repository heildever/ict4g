function x = rect_to_polar(Z)
module=sqrt(real(Z).^2+imag(Z).^2);
argument=angle(Z)*180/pi;
x=[module argument];
end