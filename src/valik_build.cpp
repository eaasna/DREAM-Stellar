#include <valik/build/build.hpp>
#include <valik/build/index_factory.hpp>
#include <valik/build/store_index.hpp>

namespace valik::app
{

void valik_build(build_arguments const & arguments)
{
    if (arguments.verbose)
    {
        std::cout << "\n-----------Index parameters-----------\n";
        std::cout << "kmer size " << std::to_string(arguments.kmer_size) << '\n';
        std::cout << "window size " << std::to_string(arguments.window_size) << '\n';
        std::cout << "bins " << std::to_string(arguments.seg_count) << '\n';
        std::cout << "bits " << std::to_string(arguments.bits) << '\n';        
        std::cout << "FPR " << std::to_string(arguments.fpr) << '\n'; 
    }

    index_factory generator{arguments};
    auto index = generator();
    store_index(arguments.out_path, index);
    return;
}

} // namespace valik::app
