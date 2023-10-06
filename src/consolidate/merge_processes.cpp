#include <utilities/consolidate/merge_processes.hpp>

namespace valik
{

bool merge_processes(search_arguments const & arguments,
                     app::search_time_statistics & time_statistics,
                     app::execution_metadata & exec_meta,
                     env_var_pack const & var_pack)
{
    // merge metadata from all threads
    exec_meta.merge(arguments, time_statistics);
    std::vector<std::string> merge_process_args;
    if (exec_meta.output_files.size() > 0)
    {
        merge_process_args.push_back(var_pack.merge_exec);
        for (auto & path : exec_meta.output_files)
            merge_process_args.push_back(path);
    }
    else
    {
        //!WORKAROUND: merge hangs if no valik matches found
        merge_process_args.push_back("echo");
        merge_process_args.push_back("-n");
    }

    external_process merge(merge_process_args);

    std::ofstream matches_out(arguments.out_file);
    matches_out << merge.cout();

    return check_external_process_success(merge_process_args, merge);
}

}   // namespace valik
