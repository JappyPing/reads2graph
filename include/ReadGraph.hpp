// ReadGraph.hpp

#ifndef READ2GRAPH_READGRAPH_HPP
#define READ2GRAPH_READGRAPH_HPP

#include "Utils.hpp"
#include "LoggingLevels.hpp"
#include "ReadWrite.hpp"
#include "MinimizerGenerator.hpp"
#include "EdgeConstructor.hpp"
#include "GraphManager.hpp"
#include "PairWiseEditDis.hpp"
#include "OMH.hpp"

#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <getopt.h>
#include <fstream>
#include <filesystem>
#include <unistd.h>
#include <seqan3/core/debug_stream.hpp> // for debug_stream

#include <array>  // std::array
#include <string> // std::string
#include <vector> // std::vector
 
#include <seqan3/alphabet/all.hpp>
#include <seqan3/alphabet/views/all.hpp> // optional: use views to convert the input string to a dna5 sequence
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/utility/views/all.hpp> // optional: use views to convert the input string to a dna5 sequence
#include <seqan3/io/sequence_file/all.hpp>

// namespace reads2graph {
    class ReadGraph{
        public:
            ReadGraph(cmd_arguments args);
            Graph reads_graph();
        private:   
            cmd_arguments args;     
    };

// }

#endif