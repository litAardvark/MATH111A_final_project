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
filename = 'data_simplified.csv';
dataTable = readtable(filename, 'ReadVariableNames', false);

%analyze_by_week(1, dataTable);
[modelOneT, modelOneD, modelOneW, modelTwoT, modelTwoD, modelTwoW] = getErrorsOverTheWeeks(dataTable);

%%
% Plot results
figure;
plot(modelOneT, 'b-', 'LineWidth', 2);
hold on;
plot(modelTwoT, 'r-', 'LineWidth', 2);
grid on;
hold off;

figure;
plot(modelOneW, 'b-', 'LineWidth', 2);
hold on;
plot(modelTwoW, 'r-', 'LineWidth', 2);
grid on;
hold off;
%%
% selectedWeek = 15;
% plotError(getErrors(getWeek(selectedWeek, dataTable), getModelTwoWeek(selectedWeek)));


% Function to get data for a given week
function result = getWeek(w, dataTable)
    if(w < 1 || w > 15)
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

% Function to compute R_out
function result = R_out(week, day, time)
    t = time-1; %time is the nth hour the Center is open
    offset = 2;%1200/600;
   result = exp(t-offset);
end


% Function to compute R_in
function result = model_two_compute(week, day, time)
    result = round(R_in(week, day, time)-R_out(week, day, time));
    if(result<0)
        result = 0;
    end
end


% Function to compute alpha
function result = getAlpha(week, day)
    %parameters
    pct = 0.08;
    AvgCsz = 26;
    AvgNumClasses = 24;
    wk8f = 1.8;
    %variables
    w = week;
    d = day; %day is the nth
    %alpha_prime = smoothAlphaPrime(w, (3*cos(d*pi)+AvgNumClasses));
    alpha_prime = 3*cos(d*pi)+AvgNumClasses;
    if (w == 8)
        pct = pct*wk8f;
    elseif (w<8)
        pct = pct*(1/5)*(log(w-0.7)+(8/exp(1)));
    elseif(w>8)
        pct = pct*(1-((1/5)*log(w-0.7)-(5/(5*exp(1)))));
    else
        pct = pct*1;
    end

    proportion = pct*AvgCsz;
    result = alpha_prime*proportion;
end

% Function to smooth alpha NOT USED
%function result = smoothAlphaPrime(week, alphaPrime)
 %   %parameters
 %   smoothedPrime = alphaPrime;
%if(week <= 4)
%    smoothedPrime = alphaPrime+((-1)^week)*((1/4)*(15+(5*(week-2.5)^2)));
%elseif (week>=8 && week<=10)
 %   smoothedPrime = alphaPrime+((-1)^(floor((week-8)/2))*(44/28));
%elseif(week>=12)
 %   smoothedPrime = alphaPrime+((1/4)*(-1)^(week-10)*(-1)^(floor((week-12)/2))*10);
%end

%result = smoothedPrime;
%end

%function to compute errors over the weeks
function [modelOne_perT, modelOne_perD, modelOne_perW, modelTwo_perT, modelTwo_perD, modelTwo_perW] = getErrorsOverTheWeeks(dataTable)
modelOne_perT = zeros(1,15);
modelOne_perD = zeros(1,15);
modelOne_perW = zeros(1,15);
modelTwo_perT = zeros(1,15);
modelTwo_perD = zeros(1,15);
modelTwo_perW = zeros(1,15);
for w=1:15
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

%Function to plot errors by week
function result = plotError(error)
    figure;
    resultRow = zeros(1,6);
    plot(resultRow);
    hold on;
    for i=1:4
        temp = error(i, :);
        for j=1:6
            resultRow(j) = temp.(j);
        end
        plot(resultRow);
    end

    hold off;
end

%Function to get model one output for given week
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

%Function to get model two output for given week
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
        error_row = array2table(wk_mod{i,:} - wk_dt{i, :});
        error_table = [error_table; error_row];
    end
    resultTable = error_table;

end

%Function to get error for given week
function [perT, perD, totalW] = getAvgError(wk_err)
    t_total_error = 0;
    d_total_error = 0;

    for i=1:4
        for j=1:6
        t_total_error = t_total_error + wk_err{i,j};
        end
        d_total_error = d_total_error + wk_err{i, 7};
    end
    perT = t_total_error/24;
    perD = d_total_error/4;
    totalW = d_total_error;
end
