#ifndef FAST_PCA_LOG_H_
#define FAST_PCA_LOG_H_

#include <cstdio>
#include <cstdlib>
#include <cstring>

#define THIS_FILE                                                       \
  (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)

#define LOG(fmt, ...)                           \
  fprintf(                                      \
      stderr, "LOG[%s:%d]: " fmt "\n",          \
      THIS_FILE, __LINE__, ##__VA_ARGS__)

#define ERROR(fmt, ...) {                                               \
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

#define CHECK_MSG(cond, fmt, ...)                                   \
  if (!(cond)) {                                                    \
    fprintf(                                                        \
        stderr, "ERROR[%s:%d]: " fmt "\n",                          \
        THIS_FILE, __LINE__, ##__VA_ARGS__);                        \
    exit(1);                                                        \
  }

#endif
