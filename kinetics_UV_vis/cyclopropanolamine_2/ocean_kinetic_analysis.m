%% load the files 
% get the files from the directory
% files = dir('20250331_mbg_5_84_660uM/*.txt');
% % convert the structured data file to table
% files = struct2table(files);
% %sort the table by column "date"
% files_sorted = sortrows(files, 'date');

%% get the data from files; write into file
% data_660 = [];
% 
% for i = 1:height(files_sorted)
%     filePath = strcat(char(files_sorted.folder(i)), '/',char(files_sorted.name(i)));
%     tempData = importdata(filePath);  % Load data
%     data_660(:,i) = tempData.data(:,2);
% end
% wavelengths = tempData.data(:,1);
% time = 1:10:7020;

data = readDataFromDir('cyclopropanolamine/20250331_mbg_5_84_660uM/*.txt');
%% plot the manually selected subsection of the wavelengths for all the spectra
% mesh(time, wavelengths(130:500), data(130:500,:))

%% get the data from files; write into file
% get the files from the directory
% files = dir('standards/*.txt');
% convert the structured data file to table
% files = struct2table(files);

[standard_spectra, wavelengths] = readDataFromDir('cyclopropanolamine_2/standards/*.txt');

% for i = 1:height(files)
%     filePath = strcat(char(files.folder(i)), '/',char(files.name(i)));
%     tempData = importdata(filePath);  % Load data
%     standards(:,i) = tempData.data(:,2);
% end

spectrum_start = 150;
trunc_standard_spectra = standard_spectra(spectrum_start:end, :);
% data = data(spectrum_start:end,:);
% wavelengths = wavelengths(spectrum_start:end,:);
start_fft_coeff = 2;
end_fft_coeff = 50;
C = [1000 800 600 400 62.5 125 250];
[M, see] = calibrationMatrix(trunc_standard_spectra, C, start_fft_coeff, end_fft_coeff);

initialRates = [];
initialConcs = [400 600 800 1000];
%% Concs vs time of 400 uM run.
[spectra, wavelengths] = readDataFromDir('cyclopropanolamine_2/20250408_mbg_5_84_400uM/*.txt');

trunc_spectra = spectra(spectrum_start:end,:);
spectra_fft = real(fft([flipud(trunc_spectra);trunc_spectra(1:end-1,:)]));
trunc_spectra_fft = spectra_fft(start_fft_coeff:end_fft_coeff,:);
concs = (M*trunc_spectra_fft).';
time = 1:10:size(spectra,2)*10;

% Compute initial rate.
start_index = 4;
concs_pts = concs(start_index:end);
time_pts = (time(start_index:end)-time(start_index)).';
stoichCoeff = 1;
rxn_percentage = 10;
[~, pts_used] = compInitialRate(stoichCoeff, time_pts, concs_pts, rxn_percentage);

A = [ones(pts_used,1) time_pts(1:pts_used)];
b = concs_pts(1:pts_used);
x = A\b;

initialRates(1) = x(2);

subplot(3,2,1)
plot(time, concs,'k.', 'LineWidth',2);
ax = gca;
ylims = ax.YLim;
hold on
plot(time, x(1)+(time - time_pts(1))*x(2), 'r', 'LineWidth',2);
ylim(ylims);
hold off
title('400 uM')
xlabel('time /s')
ylabel('concentration / uM')

%% Concs vs time of 600 uM run.
[spectra, wavelengths] = readDataFromDir('cyclopropanolamine_2/20250407_mbg_5_84_600uM/*.txt');

trunc_spectra = spectra(spectrum_start:end,:);
spectra_fft = real(fft([flipud(trunc_spectra);trunc_spectra(1:end-1,:)]));
trunc_spectra_fft = spectra_fft(start_fft_coeff:end_fft_coeff,:);
concs = (M*trunc_spectra_fft).';
time = 1:10:size(spectra,2)*10;

% Compute initial rate.
start_index = 1;
concs_pts = concs(start_index:end);
time_pts = (time(start_index:end)-time(start_index)).';
stoichCoeff = 1;
rxn_percentage = 20;
[~, pts_used] = compInitialRate(stoichCoeff, time_pts, concs_pts, rxn_percentage);

A = [ones(pts_used,1) time_pts(1:pts_used)];
b = concs_pts(1:pts_used);
x = A\b;

initialRates(2) = x(2);

subplot(3,2,2)
plot(time, concs,'k.', 'LineWidth',2);
ax = gca;
ylims = ax.YLim;
hold on
plot(time, x(1)+(time - time_pts(1))*x(2), 'r', 'LineWidth',2);
ylim(ylims);
hold off
title('600 uM')
xlabel('time /s')
ylabel('concentration / uM')

%% Concs vs time of 800 uM run.
[spectra, wavelengths] = readDataFromDir('cyclopropanolamine_2/20250407_mbg_5_84_800uM/*.txt');

trunc_spectra = spectra(spectrum_start:end,:);
spectra_fft = real(fft([flipud(trunc_spectra);trunc_spectra(1:end-1,:)]));
trunc_spectra_fft = spectra_fft(start_fft_coeff:end_fft_coeff,:);
concs = (M*trunc_spectra_fft).';
time = 1:10:size(spectra,2)*10;

% Compute initial rate.
start_index = 1;
concs_pts = concs(start_index:end);
time_pts = (time(start_index:end)-time(start_index)).';
stoichCoeff = 1;
rxn_percentage = 20;
[~, pts_used] = compInitialRate(stoichCoeff, time_pts, concs_pts, rxn_percentage);

A = [ones(pts_used,1) time_pts(1:pts_used)];
b = concs_pts(1:pts_used);
x = A\b;

initialRates(3) = x(2);

subplot(3,2,3)
plot(time, concs,'k.', 'LineWidth',2);
ax = gca;
ylims = ax.YLim;
hold on
plot(time, x(1)+(time - time_pts(1))*x(2), 'r', 'LineWidth',2);
ylim(ylims);
hold off
title('800 uM')
xlabel('time /s')
ylabel('concentration / uM')


%% Concs vs time of 1000 uM run.
[spectra, wavelengths] = readDataFromDir('cyclopropanolamine_2/20250407_mbg_5_84_1000uM/*.txt');

trunc_spectra = spectra(spectrum_start:end,:);
spectra_fft = real(fft([flipud(trunc_spectra);trunc_spectra(1:end-1,:)]));
trunc_spectra_fft = spectra_fft(start_fft_coeff:end_fft_coeff,:);
concs = (M*trunc_spectra_fft).';
time = 1:10:size(spectra,2)*10;

% Compute initial rate.
start_index = 1;
concs_pts = concs(start_index:end);
time_pts = (time(start_index:end)-time(start_index)).';
stoichCoeff = 1;
rxn_percentage = 10;
[~, pts_used] = compInitialRate(stoichCoeff, time_pts, concs_pts, rxn_percentage);

A = [ones(pts_used,1) time_pts(1:pts_used)];
b = concs_pts(1:pts_used);
x = A\b;

initialRates(4) = x(2);

subplot(3,2,4)
plot(time, concs,'k.', 'LineWidth',2);
ax = gca;
ylims = ax.YLim;
hold on
plot(time, x(1)+(time - time_pts(1))*x(2), 'r', 'LineWidth',2);
ylim(ylims);
hold off
title('1000 uM')
xlabel('time /s')
ylabel('concentration / uM')

%% Concs vs time of 334 uM run.
[spectra, wavelengths] = readDataFromDir('20250327_mbg_5_84_334uM/*.txt');

trunc_spectra = spectra(spectrum_start:end,:);
spectra_fft = real(fft([flipud(trunc_spectra);trunc_spectra(1:end-1,:)]));
trunc_spectra_fft = spectra_fft(start_fft_coeff:end_fft_coeff,:);
concs = (M*trunc_spectra_fft).';
time = 1:10:size(spectra,2)*10;

% Compute initial rate.
start_index = 2;
concs_pts = concs(start_index:end);
time_pts = (time(start_index:end)-time(start_index)).';
stoichCoeff = 1;
rxn_percentage = 20;
[~, pts_used] = compInitialRate(stoichCoeff, time_pts, concs_pts, rxn_percentage);

A = [ones(pts_used,1) time_pts(1:pts_used)];
b = concs_pts(1:pts_used);
x = A\b;

plot(time, concs,'k.', 'LineWidth',2);
ax = gca;
ylims = ax.YLim;
hold on
plot(time, x(1)+(time - time_pts(1))*x(2), 'r', 'LineWidth',2);
ylim(ylims);
hold off

%% Concs vs time of 800 uM run.
[spectra, wavelengths] = readDataFromDir('cyclopropanolamine/20250402_mbg_5_84_800uM/*.txt');

trunc_spectra = spectra(spectrum_start:end,:);
spectra_fft = real(fft([flipud(trunc_spectra);trunc_spectra(1:end-1,:)]));
trunc_spectra_fft = spectra_fft(start_fft_coeff:end_fft_coeff,:);
concs = (M*trunc_spectra_fft).';
time = 1:10:size(spectra,2)*10;

% Compute initial rate.
start_index = 1;
concs_pts = concs(start_index:end);
time_pts = (time(start_index:end)-time(start_index)).';
stoichCoeff = 1;
rxn_percentage = 10;
[~, pts_used] = compInitialRate(stoichCoeff, time_pts, concs_pts, rxn_percentage);

A = [ones(pts_used,1) time_pts(1:pts_used)];
b = concs_pts(1:pts_used);
x = A\b;

plot(time, concs,'k.', 'LineWidth',2);
ax = gca;
ylims = ax.YLim;
hold on
plot(time, x(1)+(time - time_pts(1))*x(2), 'r', 'LineWidth',2);
ylim(ylims);
hold off

%%
subplot(3,2,5:6)
plot(initialConcs, -initialRates,'ro')
title('Photochemical "rate constant"')
xlabel('Initial substrate concentration / uM')
ylabel('Initial rate / s^{-1}')
