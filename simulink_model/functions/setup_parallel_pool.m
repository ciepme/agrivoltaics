function pool = setup_parallel_pool(use_parallel)
    if nargin < 1
        use_parallel = false;
    end

    pool = [];
    if ~use_parallel
        return;
    end

    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool('local');
    end

    root_dir = pwd;
    path_def = genpath(root_dir);
    cd_futures = parfevalOnAll(pool, @cd, 0, root_dir);
    wait(cd_futures);
    path_futures = parfevalOnAll(pool, @addpath, 0, path_def);
    wait(path_futures);
end
