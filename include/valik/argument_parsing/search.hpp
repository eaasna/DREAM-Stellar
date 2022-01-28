#pragma once

#include <valik/argument_parsing/shared.hpp>

namespace valik::app
{

void init_search_parser(seqan3::argument_parser & parser, search_arguments & arguments);
void run_search(seqan3::argument_parser & parser);

} // namespace valik::app
