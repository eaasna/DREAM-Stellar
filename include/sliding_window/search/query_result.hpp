#pragma once

namespace sliding_window
{

// "query" is a noun not a verb
class query_result
{
private:
    std::string id;
    std::set<size_t> bin_hits;

public:
    query_result() = default;
    query_result(query_result const &) = default;
    query_result & operator=(query_result const &) = default;
    query_result(query_result &&) = default;
    query_result & operator=(query_result &&) = default;
    ~query_result() = default;

    query_result(std::string const query_id, std::set<size_t> query_hits)
    {
        id = query_id;
        bin_hits = query_hits; 
    }

    std::set<size_t> get_hits()
    {
        return bin_hits;
    }

    std::string get_id()
    {
        return id;
    }
    
};

} // namespace sliding_window
