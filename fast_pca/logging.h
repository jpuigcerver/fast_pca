/*
  The MIT License (MIT)

  Copyright (c) 2014,2015 Joan Puigcerver

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

#ifndef FAST_PCA_LOGGING_H_
#define FAST_PCA_LOGGING_H_

#include <cstdio>
#include <cstdlib>
#include <cstring>

#define THIS_FILE                                                       \
  (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)

#define INFO(msg)                                \
  fprintf(                                       \
      stderr, "INFO[%s:%d]: " msg "\n",          \
      THIS_FILE, __LINE__)

#define INFO_FMT(fmt, ...)                       \
  fprintf(                                       \
      stderr, "INFO[%s:%d]: " fmt "\n",          \
      THIS_FILE, __LINE__, ##__VA_ARGS__)

#define WARN(msg)                                \
  fprintf(                                       \
      stderr, "WARN[%s:%d]: " msg "\n",          \
      THIS_FILE, __LINE__)

#define WARN_FMT(fmt, ...)                       \
  fprintf(                                       \
      stderr, "WARN[%s:%d]: " fmt "\n",          \
      THIS_FILE, __LINE__, ##__VA_ARGS__)


#define ERROR(msg) {                                                    \
    fprintf(                                                            \
        stderr, "ERROR[%s:%d]: " msg "\n",                              \
        THIS_FILE, __LINE__);                                           \
    exit(1);                                                            \
  }

#define ERROR_FMT(fmt, ...) {                                           \
    fprintf(                                                            \
        stderr, "ERROR[%s:%d]: " fmt "\n",                              \
        THIS_FILE, __LINE__, ##__VA_ARGS__);                            \
    exit(1);                                                            \
  }

#define CHECK(cond)                                                     \
  if (!(cond)) {                                                        \
    fprintf(                                                            \
        stderr, "ERROR[%s:%d]: Check failed (" #cond ")\n",             \
        THIS_FILE, __LINE__);                                           \
    exit(1);                                                            \
  }

#define CHECK_MSG(cond, msg)                                     \
  if (!(cond)) {                                                 \
    fprintf(                                                     \
        stderr, "ERROR[%s:%d]: " msg "\n", THIS_FILE, __LINE__); \
    exit(1);                                                     \
  }

#define CHECK_FMT(cond, fmt, ...)                                   \
  if (!(cond)) {                                                    \
    fprintf(                                                        \
        stderr, "ERROR[%s:%d]: " fmt "\n",                          \
        THIS_FILE, __LINE__, ##__VA_ARGS__);                        \
    exit(1);                                                        \
  }

#endif  // FAST_PCA_LOGGING_H_
