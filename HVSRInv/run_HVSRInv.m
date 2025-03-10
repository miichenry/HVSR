% Set up the initial parameters
LayerNum = 2;               % Number of layers
InvNum = 5;                 % Number of inversions
IterNum = 25000;             % Number of iterations
PopNum = 1000;               % Population size
GeneNum = 20;               % Number of genes
% Set additional parameters
FreqMin = 0.3;              % Minimum frequency for analysis
FreqMax = 48;              % Maximum frequency for analysis
%SampleNum = 1000;            % Number of samplesi

% Load frequency and HVSR data from a text file (first and second columns)
% Check if filename variable exists
if ~exist('filename', 'var')
    error('Please provide the filename as an argument.');
end

% Display the filename
disp('Initial Parameters:');
disp(repmat('-', 1, 120));
disp(['LayerNum: ', num2str(LayerNum)]);
disp(['InvNum: ', num2str(InvNum)]);
disp(['IterNum: ', num2str(IterNum)]);
disp(['PopNum: ', num2str(PopNum)]);
disp(['GeneNum: ', num2str(GeneNum)]);
disp(['FreqMin: ', num2str(FreqMin)]);
disp(['FreqMax: ', num2str(FreqMax)]);
disp(['Using data from file: ', filename]);
data = readmatrix(filename); % Read the entire data from the file

% Extract number from the filename using regular expression
match = regexp(filename, '(\d+)_.*\.txt', 'once', 'tokens');

% If a match is found, extract the number
if ~isempty(match)
    sta_name = match{1};
    disp(['Station: ', sta_name]);
else
    disp('No station name found in filename.');
end

% Extract the first column for frequency and second column for HVSRData
Freq = data(:, 1);          % Extract the first column for frequency
HVSRData = data(:, 2);      % Extract the second column for HVSRData

% Check that HVSR data is not empty
if isempty(Freq) || isempty(HVSRData)
    error('HVSR data is empty. Please check the data file.');
end

format long g;
% Load model parameters (InitModData) from another text file (format as specified)
model_params_filename = 'subsurface1.txt'; % Replace with your model parameters file name
model_data = readmatrix(model_params_filename);  % Read all data from the model parameters file

% Open the file for reading
fileID = fopen(model_params_filename, 'r');

% Initialize InitModData array
InitModData = zeros(LayerNum, 8);  % Initialize the matrix for model parameters

% Read the file line by line
layer_index = 1;  % To keep track of layers

while ~feof(fileID)
    % Read the Layer header line
    layer_header = fgetl(fileID);
    
    % Read the parameters for the current layer
    if layer_header ~= -1
        param_line = fgetl(fileID);  % Read the next line which contains the parameters
        params = str2double(strsplit(param_line));  % Convert the string to numbers
        
        % Store the parameters into InitModData
        if length(params) == 8
            InitModData(layer_index, :) = params;  % Store the 8 parameters for the current layer
            layer_index = layer_index + 1;  % Move to the next layer
        else
            warning('Expected 8 parameters for layer %d, but got %d. Skipping layer.', layer_index, length(params));
        end
    end
end

% Close the file after reading
fclose(fileID);

% Check InitModData to ensure it is populated correctly
if isempty(InitModData)
    error('InitModData is empty. Please check the model parameter file.');
end

% Display the loaded InitModData for confirmation
disp('Loaded InitModData:');
disp(InitModData);

[minValue,closestIndex1] = min(abs(Freq-FreqMin));
[minValue,closestIndex2] = min(abs(Freq-FreqMax));
SampleNum=closestIndex2-closestIndex1+1;
disp(['FreqMin:FreqMax:  ', num2str(FreqMin), ':', num2str(FreqMax)]);
disp(['SampleNum: ', num2str(SampleNum), ', closestIndex2: ', num2str(closestIndex2), ', closestIndex1: ', num2str(closestIndex1)]);

% Check FreqMin and FreqMax values
if FreqMin <= 0 || FreqMax <= 0
    error('Frequency values must be greater than 0.');
end

if FreqMin >= FreqMax
    error('FreqMin must be less than FreqMax.');
end

% Check LayerNum value
if LayerNum < 2
    error('LayerNum must be at least 2.');
end

% Check that frequency range is within the provided data range
if Freq(1) > FreqMin
    error('The minimum frequency can not be less than %g', Freq(1));
end

if Freq(end) < FreqMax
    error('The maximum frequency can not be greater than %g', Freq(end));
end

% Create an "Outputs" directory if it doesn't exist
if ~exist('Outputs', 'dir')
    mkdir('Outputs');
end

% Create output.txt file for saving inversion details
[~, baseFileName, ~] = fileparts(filename);
outputFileBase = baseFileName;
outputFile = fullfile('Outputs', [sta_name '_output.txt']);
fileID = fopen(outputFile, 'wt');
fprintf(fileID, 'HVSR data file: %s\n', filename);
fprintf(fileID, 'Inversion Parameters\n');
fprintf(fileID, 'Number_of_Layers: %d\n', LayerNum);
fprintf(fileID, 'Minimum_Frequency: %2.3f\n', FreqMin);
fprintf(fileID, 'Maximum_Frequency: %2.3f\n', FreqMax);
fprintf(fileID, 'Number_of_Genes: %d\n', GeneNum);
fprintf(fileID, 'Number_of_Populations: %d\n', PopNum);
fprintf(fileID, 'Number_of_Iterations: %d\n', IterNum);
fprintf(fileID, 'Number_of_Inversions: %d\n', InvNum);

% Write initial model data to output file
for i = 1:LayerNum
    S = InitModData(i, :);
    fprintf(fileID, 'Layer_No: %d\n', i);
    fprintf(fileID, 'HMin HMax VMin VMax DenMin DenMax DampMin DampMax\n');
    fprintf(fileID, '%1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f\n', S(1), S(2), S(3), S(4), S(5), S(6), S(7), S(8));
end

% Start timer and run Inversion process
beginpar = tic;
par = Inversion(sta_name, LayerNum, InvNum, IterNum, PopNum, GeneNum, Freq, HVSRData, FreqMin, FreqMax, SampleNum, InitModData, gca);
endpar = toc(beginpar);
fprintf('Inversion took: %f seconds\n', endpar);

% Get the size of the inversion results (par)
sz = size(par);

% Write inversion results to output file if direc == 1        
for i = 1:sz(1)
    fprintf(fileID, 'Inversion_No:%d\n', i);
    S = par(i, :);   
    fprintf(fileID, 'H V Den Damp\n');
    a = 1;
    b = 4;
    for j = 1:sz(2) / 4   
        c = S(a:b);
        fprintf(fileID, '%1.4f %1.4f %1.4f %1.4f\n', c(1), c(2), c(3), c(4));
        a = a + 4;
        b = b + 4;
    end
    fprintf(fileID, 'Goodness of Fit: %4f\n', S(a));  % Assuming goodness of fit is stored after the parameters
end;

for i=1:sz(2)
        hMean(1,i)=mean(par(:,i));
end;

% Write the average model to the output file if direc == 1
fprintf(fileID, 'Average Model\n');
fprintf(fileID, 'H V Den Damp\n');

% Write the averaged parameters (H, Vs, Den, Damp) to the output file
j=1;i=1;
while i<sz(2)
    H(j)=round(hMean(i));i=i+1;
    Vs(j)=round(hMean(i));i=i+1;
    Den(j)=round(hMean(i),4);i=i+1;
    Damp(j)=round(hMean(i),4);i=i+1;      
    fprintf(fileID,'%1.4f %1.4f %1.4f %1.4f\n',H(j),Vs(j),Den(j),Damp(j));
    j=j+1;
end;
  
fclose(fileID);

% Calculate depth range for the model plot
depthmax = sum(InitModData(:, 2));

% Generate the model plot
[modelx, modely] = SetArray(H, Vs, depthmax);

% Plot the model
fig2 = figure;
plot(modelx, modely, 'LineWidth', 4, 'Color', 'r');
set(gca, 'YDir', 'reverse');
xlabel('Velocity (ms^{-1})');
ylabel('Depth (m)');
xlim([0 max(modelx)]);
ylim([0 modely(end)]);
set(gca, 'FontSize', 14);
grid on;

% Automatically save the model plot
saveas(fig2, fullfile('Outputs', [sta_name '_Vs.png']));
close(fig2);
% 
% for i=1:layerNum
%     if i<layerNum
%         results(i,1)=H(i);
%     else
%         results(i,1)=0;
%     end;
%     results(i,2)=Vs(i);
%     results(i,3)=Den(i);
%     results(i,4)=Damp(i);    
% end;

% Calculate the estimated HVSR using the model
frecInt = FreqMax / SampleNum;  % Frequency interval for the new frequency array
f1 = FreqMin;
f2 = FreqMax;
freqs2 = (f1:frecInt:f2);  % New frequency range for the estimated HVSR

% Calculate the estimated HVSR using the parameters (model output)
HVSR2 = CalcHVSR(Vs', H', Den', Damp', freqs2);  % Assuming Vs, H, Den, Damp are defined

% Plot both original and estimated HVSR
fig1 = figure;
semilogx(Freq, HVSRData, 'k', 'LineWidth', 1.5);  % Original HVSR
hold on;
semilogx(freqs2, HVSR2, '--r', 'LineWidth', 2);   % Estimated HVSR
hold off;

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('Amplitude');
xlim([FreqMin,FreqMax]);
set(gca, 'FontSize', 14);
legend('HVSR', 'Estimated HVSR', 'Location', 'northwest');
grid on;

% Automatically save the HVSR plot
saveas(fig1, fullfile('Outputs', [sta_name '_HVSR.png']));
close(fig1);
