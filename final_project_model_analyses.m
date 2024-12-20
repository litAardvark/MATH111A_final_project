% MATH 111A Final Project
% Author: Kwadwo Danquah
% Date: 12NOV24
% Description: This script tests my model against the underlying data and
% plots errors

%% Initialization
clear;          % Clear variables
clc;            % Clear command window
close all;      % Close all figure windows

%% Main Script
% Example inputs
filename = 'data_simplified.csv';
dataTable = readtable(filename, 'ReadVariableNames', false);

% % % Call functions
[modelOneT, modelOneD, modelOneW, modelTwoT, modelTwoD, modelTwoW] = getErrorsOverTheWeeks(dataTable);
% 
% %%
% % Plot results
figure;
% plot(modelOneT, 'b-', 'LineWidth', 2);
% hold on;
colors = getColor(modelTwoT, 5); 
b = bar(modelTwoT);
b.FaceColor = 'flat';
b.CData = colors;
hold off;
title('mae per hour');
xlabel('week');
ylabel('error/t');

figure;
% plot(modelOneD, 'b-', 'LineWidth', 2);
% hold on;
colors = getColor(modelTwoD, 10); 
b = bar(modelTwoD);
b.FaceColor = 'flat';
b.CData = colors;
hold off;
title('error per daily count');
xlabel('week');
ylabel('error/d');
%
% for i=1:16
%  selectedWeek = i;
%  plotOutput(getModelTwoWeek(selectedWeek), getWeek(selectedWeek, dataTable), selectedWeek);
%  plotError(getErrors(getWeek(selectedWeek, dataTable), getModelTwoWeek(selectedWeek)), selectedWeek);
% end


%% Function Definitions
%Function to get Bar Colors
function result = getColor(data, threshold)
colors = zeros(length(data), 3);  % Initialize color matrix
for i = 1:length(data)
    if abs(round(data(i))) > threshold
        colors(i, :) = [1, 0, 0];  % Red for values above threshold
    else
        colors(i, :) = [0, 0, 1];  % Blue for values below threshold
    end
end
result = colors;
end





% Function to get data for a given week
function result = getWeek(w, dataTable)
    if(w < 1 || w > 16)
        w = 1;
    end
    result_table = table();
    numRows = 4;
    curSect = ((w-1)*numRows);
    for i = 1:numRows
        sourceRow = dataTable((curSect+i), :);
        result_table = [result_table; sourceRow]; %appending row from source
    end
    result = result_table;
end

% Function to compute model 1 output
function [out1, out2] = model_one_compute(week, day, time)
    divisor = 6;
    t = time-1; %time is the nth hour the Center is open
    alpha_prime = getAlpha(week, day);
    if(t<2)
        timeCount = round((alpha_prime/divisor)*((-2/pi)*cos((t*pi)/2)+t+(pi+2)/pi));
    else
        timeCount = round((alpha_prime/divisor)*((2/pi)*(cos((t-2)*pi/2)-cos((t*pi)/2))+4));
    end
    total = round(alpha_prime);


    out1 = timeCount;
    out2 = total;
end

% Function to compute R_in
function result = R_in(week, day, time)
    t = time-1; %time is the nth hour the Center is open
    divisor = 6;
    alpha_prime = getAlpha(week, day);
    timeCount = (alpha_prime/divisor)*((-2/pi)*cos((t*pi)/2)+t+((pi+2)/pi));
    
    result = timeCount;
end


% Function to compute R_in
function result = R_in_piecewise(week, day, time)
    t = time-1; %time is the nth hour the Center is open
    [alpha_prime, beta] = getAlpha(week, day);
    if t<3
        result = (alpha_prime/6)*((-1/pi)*cos(t*pi/2)+t+((pi+1)/pi));
    else
        %result = (0.633.^(t-3))*(alpha_prime/24)*((-10/pi)*cos((t-3)*pi/2)+(2*(t-3))+((16*pi+14)/pi));
        result = (alpha_prime/24)*((-10/pi)*cos((t-3)*pi/2)+(2*(t-3))+((16*pi+14)/pi));
    end
    
end



% Function to compute R_out
function result = R_out(week, day, time)
    t = time-1; %time is the nth hour the Center is open
    offset = 2;%1200/600;
   result = exp(t-offset);
end

% Function to compute R_out Technical
function result = R_out_technical(week, day, time)
    t = time-1; %time is the nth hour the Center is open
    result = 0;
   for n = 0:(t-1)
       R_i = R_in_piecewise(week,day,(1));
       result = result + R_i*exp(-(t/0.5));
   end
end


% Function to compute R_in
function result = model_two_compute(week, day, time)
    %result = round(R_in(week, day, time)-R_out(week, day, time));
    result = round(R_in_piecewise(week, day, time)-R_out_technical(week, day, time));
    if(result<0)
        result = 0;
    end
end


% Function to compute alpha
function [alpha_prime, beta] = getAlpha(week, day)
    %parameters
    pct = 0.08;
    AvgCsz = 26;
    AvgNumClasses = 24;
    wk8f = 1.8;
    %variables
    w = week;
    d = day; %day is the nth
    %alpha_prime = smoothAlphaPrime(w, (3*cos(d*pi)+AvgNumClasses));
    beta = 3*cos(d*pi)+AvgNumClasses;
    if (w == 8 || w==16)
        pct = pct*wk8f;
    elseif (w<8)
        pct = pct*(1/5)*(log(w-0.7)+(8/exp(1)));
    elseif(w>8)
        pct = pct;%pct*(1-((1/5)*log(w-0.7)-(5/(5*exp(1)))));
    else
        pct = pct*1;
    end

    proportion = pct*AvgCsz;
    alpha_prime = beta*proportion;
end


%function to compute errors over the weeks
function [modelOne_perT, modelOne_perD, modelOne_perW, modelTwo_perT, modelTwo_perD, modelTwo_perW] = getErrorsOverTheWeeks(dataTable)
modelOne_perT = zeros(1,16);
modelOne_perD = zeros(1,16);
modelOne_perW = zeros(1,16);
modelTwo_perT = zeros(1,16);
modelTwo_perD = zeros(1,16);
modelTwo_perW = zeros(1,16);
for w=1:16
    data = getWeek(w, dataTable);
    disp(['WEEK ', num2str(w)]);
    disp('data:');
    disp(data);
    model = getModelOneWeek(w);
    error = getErrors(data, model);
    [avgPerT, avgPerD, perW] = getAvgError(error);
    modelOne_perT(w) = avgPerT;
    modelOne_perD(w) = avgPerD;
    modelOne_perW(w) = perW;
    disp('model one:');
    disp(model);
    disp('model one error:');
    disp(error);
    model = getModelTwoWeek(w);
    error = getErrors(data, model);
    [avgPerT, avgPerD, perW] = getAvgError(error);
     modelTwo_perT(w) = avgPerT;
    modelTwo_perD(w) = avgPerD;
    modelTwo_perW(w) = perW;
    disp('model two:');
    disp(model);
    disp('model two error:');
    disp(error);
end
end


% customColors = [0 0.4470 0.7410;    % Blue
%                 0.8500 0.3250 0.0980; % Red
%                 0.9290 0.6940 0.1250; % Yellow
%                 0.4940 0.1840 0.5560; % Purple
%                 0.4660 0.6740 0.1880; % Green
%                 0.3010 0.7450 0.9330; % Cyan
%                 0.6350 0.0780 0.1840]; % Dark Red
% set(groot, 'defaultAxesColorOrder', customColors);
% repeatedColors = [customColors; customColors];

%Function to plot errors by week
function result = plotError(error, selectedWeek)
    figure;
    
    resultRow = zeros(1,6);
    plot(resultRow);
    hold on;
    for i=1:4
        temp = error(i, :);
        for j=1:6
            resultRow(j) = temp.(j);
        end
        plot(resultRow, 'LineWidth', 1.5);
    end

    hold off;
    ax = gca;
    ax.XTickLabel = {'0', '1', '2', '3', '4', '5'};
    titleString = ['week', ' ', num2str(selectedWeek), ' ', 'error'];
    title(titleString);
    legend('', 'Mon', 'Tue', 'Wed', 'Thu', 'Location', 'best');
    xlabel('t');
    ylabel('error');

end

%Function to plot outpu by week
function result = plotOutput(model, data, selectedWeek)
    modelRow = zeros(1,6);
    dataRow = zeros(1,6);
    figure; 
    customColors = [0 0.4470 0.7410;    % Blue
                0.8500 0.3250 0.0980; % Red
                0.9290 0.6940 0.1250; % Yellow
                0.4940 0.1840 0.5560; % Purple
                0.4660 0.6740 0.1880; % Green
                0.3010 0.7450 0.9330; % Cyan
                0.6350 0.0780 0.1840]; % Dark Red
    repeatedColors = [customColors; customColors];
    ax = gca;
    plot(modelRow);
    hold on;

    for i=1:4
        tempMod = model(i, :);
        tempDat = data(i,:);
        for j=1:6
            modelRow(j) = tempMod.(j);
            dataRow(j) = tempDat.(j);
        end
        ci = i+1;
        plot(modelRow, 'Color', repeatedColors(ci, : ), 'LineWidth', 2);
        plot(dataRow, 'Color', repeatedColors(ci+7, : ), 'LineWidth', 2, 'LineStyle',':');
    end
    hold off;
    ax = gca;
    ax.XTickLabel = {'0', '1', '2', '3', '4', '5'};
    titleString = ['week', ' ', num2str(selectedWeek), ' ', 'model vs data'];
    title(titleString);
    legend('', 'Mon (model)', 'Mon (data)', 'Tue (model)', 'Tue (data)', 'Wed (model)', 'Wed (data)', 'Thu (model)', 'Thu (data)', 'Location', 'best');
    xlabel('t');
    ylabel('count');

end

%Function to get model output for given week
function result = getModelOneWeek(week)
    result_table = table();
    numRows = 4;
    for d = 1:4
        row = zeros(1, 7);
        total = 0;
        for i = 1:6
            t = i;
            [output, model_total] = model_one_compute(week, d, t);
            row(i) = output;
            total = model_total;
        end
        row(7) = total;
        rowTable = array2table(row);
        result_table = [result_table;rowTable];
    end

    result = result_table;
end

%Function to get model output for given week
function result = getModelTwoWeek(week)
    result_table = table();
    numRows = 4;
    for d = 1:4
        row = zeros(1, 7);
        for i = 1:6
            t = i;
            output = model_two_compute(week, d, t);
            row(i) = output;
        end
        row(7) = round(getAlpha(week, d));
        rowTable = array2table(row);
        result_table = [result_table;rowTable];
    end

    result = result_table;
end

%Function to get error for given week
function resultTable = getErrors(wk_dt, wk_mod)
    error_table = table();

    for i=1:4
        error_row = (array2table(wk_mod{i,:} - wk_dt{i, :}));
        error_table = [error_table; error_row];
    end
    resultTable = error_table;

end

%Function to get errors for given week
function [perT, perD, totalW] = getAvgError(wk_err)
    t_total_error = 0;
    d_total_error = 0;

    for i=1:4
        for j=1:6
        t_total_error = t_total_error + abs(wk_err{i,j});
        end
        d_total_error = d_total_error + abs(wk_err{i, 7});
    end
    perT = t_total_error/24;
    perD = d_total_error/4;
    totalW = d_total_error;
end
