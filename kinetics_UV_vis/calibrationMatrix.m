function [M, see] = calibrationMatrix(spectra, concs, start_fft_coeff, end_fft_coeff)
    spectra_fft = real(fft([flipud(spectra);spectra(1:end-1,:)]));
    standards_fft_trunc = spectra_fft(start_fft_coeff:end_fft_coeff,:);
    
    A = standards_fft_trunc*(standards_fft_trunc.');
    [S,E] = eig(A);
    S = fliplr(S); e = flipud(diag(E));
    s = zeros(size(e));
    for i = 1:size(s,1)
        s(i) = e(i)/sum(e(i:end));
    end
    [~,i] = min(s);
    plot(s);
    num_eigenvecs = i-1;
    V = S(:,1:num_eigenvecs);
    Z = V.'*standards_fft_trunc;
    P = (Z.'\concs.').';
    M = P*V.';
    see = norm(concs - M*standards_fft_trunc);
end

