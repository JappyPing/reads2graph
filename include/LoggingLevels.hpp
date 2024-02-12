/**
 * @Author: Pengyao Ping
 * @Date:   2023-02-09 10:06:30
 * @Last Modified by:   Pengyao Ping
 * @Last Modified time: 2023-02-09 23:46:49
 */
#ifndef LOGGINGLEVELS_HPP
#define LOGGINGLEVELS_HPP

#include <iostream>
#include <string>

#define LOG_LEVEL_DEBUG 1
#define LOG_LEVEL_INFO 2
#define LOG_LEVEL_WARNING 3
#define LOG_LEVEL_ERROR 4
#define LOG_LEVEL_CRITICAL 5

// Set the logging level globally
extern int LOG_LEVEL;

// Define colors for different logging levels
const std::string COLOR_DEBUG = "\033[32m";   // green
const std::string COLOR_INFO = "\033[36m";    // cyan
const std::string COLOR_WARNING = "\033[33m"; // yellow
const std::string COLOR_ERROR = "\033[31m";   // red
const std::string COLOR_CRITICAL = "\033[35m";// purple
const std::string COLOR_RESET = "\033[0m";

#endif