%---------------------------- Huffman Algorithm --------------------------%
clear;
clc;
close all;

%--------------------------- Input Section -------------------------------%
% Ask user to enter probabilities for which to generate Huffman code
% Our Example: [0.35 0.30 0.20 0.10 0.04 0.005 0.005]
% Note: The sum of probabilities must equal 1
Input_Message = input('Please Enter The Given Probabilities: ');

%------------------------- Initialize Arrays -----------------------------%
% Create a zero matrix 'Arr' with:
%   rows = number of symbols
%   columns = number of symbols - 1
% Store the input probabilities (transposed) in the first column of Arr.
% Store the length of the input vector in 'Length'.
Arr = zeros(length(Input_Message), length(Input_Message) - 1);
Arr(:, 1) = transpose(Input_Message);
Length = length(Input_Message);

%-------------------------- Build Probability Tree -----------------------%
% Arrange probabilities in descending order and iteratively combine the
% smallest two until one root node remains.
for i = 1:length(Input_Message)
    Arr(:, i) = sort(Arr(:, i), 'descend');
    
    % Stop when reaching the last column
    if i == length(Input_Message) - 1
        break;
    end

    % Combine smallest two probabilities
    Arr(1:Length - 2, i + 1) = Arr(1:Length - 2, i);
    Arr(Length - 1, i + 1) = Arr(Length, i) + Arr(Length - 1, i);
    Length = Length - 1;
end

%----------------------- Initialize Huffman Code Table -------------------%
% Determine the size of Arr
N = size(Arr);
P = N(1, 2);    % Number of columns

% Initialize counters
L = 3;
M = 2;

% Create a string array for Huffman codes
Huffman_Source_Code = strings(size(Arr));

% Initialize last column codes ("0" and "1")
Huffman_Source_Code(1, P) = "0";
Huffman_Source_Code(2, P) = "1";

%------------------------ Huffman Encoding Loop --------------------------%
% Build Huffman codes column by column from right to left.
for k = P - 1:-1:1
    q = -1;
    
    % Search for matching combined probabilities
    for u = 1:M
        if (Arr(L, k) + Arr(L - 1, k) == Arr(u, k + 1))
            q = u;
            Huffman_Source_Code(L - 1, k) = strcat(Huffman_Source_Code(u, k + 1), "0");
            Huffman_Source_Code(L, k)     = strcat(Huffman_Source_Code(u, k + 1), "1");
            break;
        end
    end

    % Update codes for the remaining symbols
    j = 1;
    for y = 1:L - 2
        if (y == q)
            j = j + 1;
        end
        Huffman_Source_Code(y, k) = Huffman_Source_Code(j, k + 1);
        j = j + 1;
    end
    
    % Increment counters
    M = M + 1;
    L = L + 1;
end

%--------------------------- Display Results -----------------------------%
disp('Sorted Probabilities (Arr):');
disp(Arr);
disp('Generated Huffman Codes:');
disp(Huffman_Source_Code);

%------------------ Entropy, Average Length, Efficiency ------------------%
% Entropy = Σ P(i) * log2(1/P(i))
Entropy_of_symbols = 0;
for i = 1:length(Input_Message)
    Entropy_of_symbols = Entropy_of_symbols + (-Input_Message(i) * log2(Input_Message(i)));
end
disp(['Entropy: H(X) = ', num2str(Entropy_of_symbols), ' bits/symbol']);

% Average code length = Σ P(i) * L(i)
Avg_Len = 0;
for i = 1:length(Input_Message)
    Avg_Len = Avg_Len + length(Huffman_Source_Code{i, 1}) * Input_Message(i);
end
disp(['Average length: L(X) = ', num2str(Avg_Len), ' bits/symbol']);

% Efficiency = (Entropy / Average length) * 100%
efficiency = (Entropy_of_symbols / Avg_Len) * 100;
disp(['Efficiency = ', num2str(efficiency), '%']);

%*************************************************************************%
