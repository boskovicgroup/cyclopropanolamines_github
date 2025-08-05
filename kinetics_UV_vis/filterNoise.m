function filtered = filterNoise(intensities, num_coeffs)
    n = size(intensities,1);
    ints_fft = real(fft([flipud(intensities);intensities(1:end-1)]));
    ints_ifft = ifft([ints_fft(1:num_coeffs);zeros(2*(n-num_coeffs)-1,1);ints_fft(end-(num_coeffs-1):end)],2*n-1);
    filtered = flipud(real(ints_ifft(1:n)));
end