clc;
clear;
close all;

%%
%n = input('Enter number of symbols: ');
%symbols = cell(1, n);
%for i = 1:n
%    symbols{i} = input(sprintf('Enter symbol %d: ', i), 's');
%end

%p = zeros(1, n);
%for i = 1:n
%    p(i) = input(sprintf('Enter probability of %s: ', symbols{i}));
%end
%%%
%% Our Example
symbols = {'a', 'b', 'c', 'd', 'e', 'f', 'g'}; 
p   = [0.35, 0.3, 0.2, 0.1, 0.04, 0.005, 0.005]; 


p = p / sum(p);
[p, idx] = sort(p, 'descend');
symbols = symbols(idx);

fano_encode = @(sym, prob) recursive_fano(sym, prob);
codes = fano_encode(symbols, p);

fprintf('\nSymbol   Probability   Fano Code\n');
for i = 1:length(symbols)
    fprintf('  %-6s    %.3f        %s\n', symbols{i}, p(i), codes{i});
end

H = -sum(p .* log2(p));
L = 0;
for i = 1:length(codes)
    L = L + p(i) * length(codes{i});
end
eta = H / L;

fprintf('\nEntropy (H): %.4f bits/symbol\n', H);
fprintf('Average Code Length (L): %.4f bits/symbol\n', L);
fprintf('Coding Efficiency (η = H/L): %.2f%%\n', eta * 100);


%% ******************** Fano Code Generation Function ******************
function codes = recursive_fano(symbols, p)
    n = length(symbols);
    codes = cell(1, n);
    if n == 1
        codes{1} = '';
        return;
    end
    total = sum(p);
    diff = abs(cumsum(p) - total/2);
    [~, split] = min(diff);
    left = recursive_fano(symbols(1:split), p(1:split));
    right = recursive_fano(symbols(split+1:end), p(split+1:end));
    for i = 1:length(left)
        codes{i} = ['0' left{i}];
    end
    for i = 1:length(right)
        codes{split+i} = ['1' right{i}];
    end
end