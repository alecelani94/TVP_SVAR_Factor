function [Data, dates, names, tcode] = prepare_data(xlsx_file)
% prepare_data  Load and transform the 14-variable monthly dataset.
%
%   [Data, dates, names, tcode] = prepare_data(xlsx_file)
%
%   Replicates the data-preparation steps performed inside realVARDGP.m of
%   Korobilis (2022)'s replication package, but exposes them as a stand-alone
%   utility so the resulting matrix Data can be fed directly into any
%   econometric routine (OLS VAR, Gibbs sampler, HMC, ...).
%
%   Steps:
%     1. Read raw monthly levels (B2:O518 by default), the tcode row, and
%        the date column from the 'data' sheet of the xlsx.
%     2. Apply transx() variable-by-variable using the tcode in row 2.
%     3. Outlier-adjust each transformed series via adjout(., 4.5, 4).
%     4. Drop the first observation (lost to log-differencing).
%     5. Standardize (zscore) every column.
%
%   Inputs
%     xlsx_file : path to data_for_MC.xlsx (default: this folder).
%
%   Outputs
%     Y     : T-by-14 matrix of standardized, stationary series.
%     dates : T-by-1 datetime vector aligned with Y.
%     names : 1-by-14 cellstr of variable names.
%     tcode : 1-by-14 integer transformation codes.
%
%   The function depends on transx.m and adjout.m, which live in
%   ../Replication (SVAR Factor)/MONTE_CARLO/functions/. That folder is
%   added to the path automatically if not already on it.

    if nargin < 1 || isempty(xlsx_file)
        here      = fileparts(mfilename('fullpath'));
        xlsx_file = fullfile(here, 'data_for_MC.xlsx');
    end

    % Make Korobilis's helper functions visible
    here    = fileparts(mfilename('fullpath'));
    funcdir = fullfile(here, '..', 'Replication (SVAR Factor)', ...
                       'MONTE_CARLO', 'functions');
    if exist(funcdir, 'dir') && ~contains(path, funcdir)
        addpath(funcdir);
    end

    %% 1. Read raw block, header, dates
    hdr   = readcell(xlsx_file, 'Sheet', 'data', 'Range', 'B1:O1');
    names = hdr(1, :);                                                   % 1 x 14 cellstr
    tcode = readmatrix(xlsx_file, 'Sheet', 'data', 'Range', 'B2:O2');    % 1 x 14
    data  = readmatrix(xlsx_file, 'Sheet', 'data', 'Range', 'B3:O518');  % 516 x 14
    date_serials = readmatrix(xlsx_file, 'Sheet', 'data', 'Range', 'A3:A518');
    dates_full   = datetime(date_serials, 'ConvertFrom', 'excel');

    %% 2. Apply transformations variable-by-variable
    [T_raw, n] = size(data);
    Dt = nan(T_raw, n);
    for i = 1:n
        Dt(:, i) = transx(data(:, i), tcode(i));
        Dt(:, i) = adjout(Dt(:, i), 4.5, 4);
    end

    %% 3. Drop first row (lost to differencing) and standardize
    Data  = zscore(Dt(2:end, :));
    dates = dates_full(2:end);
end
