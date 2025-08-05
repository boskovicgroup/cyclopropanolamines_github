function [spectra, wavelengths] = readDataFromDir(path)
    % get the files from the directory
    files = dir(path);
    % convert the structured data file to table
    files = struct2table(files);
    %sort the table by column "date"
    files_sorted = sortrows(files, 'date');
    disp(files_sorted);
    spectra = [];

    for i = 1:height(files_sorted)
        filePath = strcat(char(files_sorted.folder(i)), '/',char(files_sorted.name(i)));
        tempData = importdata(filePath);  % Load data
        spectra = cat(2, spectra, tempData.data(:,2));
    end
    wavelengths = tempData.data(:,1);
end

