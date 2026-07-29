// This file is part of TetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2018 Jeremie Dumas <jeremie.dumas@ens-lyon.org>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Created by Jeremie Dumas on 09/04/18.
//

#pragma once

#include <fstream>
#include <mutex>
#include <ostream>
#include <sstream>
#include <string>

namespace floatTetWild {

// The slice of spdlog this project used: a level filter, a console sink, an optional file sink, and
// {} substitution. Every format string in the tree is a bare {} with no format spec and every
// argument streams, so writing one argument per {} covers all of them.
class Logger
{
  public:
    // Same numbering spdlog used, so --level keeps its meaning.
    enum Level
    {
        Trace = 0,
        Debug,
        Info,
        Warn,
        Error,
        Critical,
        Off
    };

    // By default, write to stdout, but don't write to any file
    static void init(bool               use_cout = true,
                     const std::string& filename = std::string(),
                     bool               truncate = true);

    void set_level(int level);

    template<typename... Args>
    void debug(const char* fmt, const Args&... args)
    { log(Debug, fmt, args...); }
    template<typename... Args>
    void info(const char* fmt, const Args&... args)
    { log(Info, fmt, args...); }
    template<typename... Args>
    void warn(const char* fmt, const Args&... args)
    { log(Warn, fmt, args...); }
    template<typename... Args>
    void error(const char* fmt, const Args&... args)
    { log(Error, fmt, args...); }

  private:
    template<typename... Args>
    void log(Level level, const char* fmt, const Args&... args)
    {
        if (level < level_)
            return;
        std::ostringstream message;
        format(message, fmt, args...);
        write(level, message.str());
    }

    // Copies fmt through, writing the next argument at each {}. Anything else, a lone brace
    // included, goes through untouched. Arguments left over once the placeholders run out are
    // dropped, where fmt would have thrown.
    static void format(std::ostream& out, const char* fmt) { out << fmt; }

    template<typename T, typename... Rest>
    static void format(std::ostream& out, const char* fmt, const T& value, const Rest&... rest)
    {
        for (const char* c = fmt; *c != '\0'; ++c) {
            if (c[0] == '{' && c[1] == '}') {
                out << value;
                format(out, c + 2, rest...);
                return;
            }
            out << *c;
        }
    }

    void write(Level level, const std::string& message);

    int           level_    = Info;  // spdlog's default
    bool          use_cout_ = true;
    std::ofstream file_;
    std::mutex    mutex_;
};

// Retrieve current logger, or create one if not available
Logger& logger();

}  // namespace floatTetWild
