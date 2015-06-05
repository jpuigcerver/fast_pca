/*
  The MIT License (MIT)

  Copyright (c) 2015 Joan Puigcerver

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

#ifndef FAST_PCA_ENDIAN_H_
#define FAST_PCA_ENDIAN_H_

#if defined(__linux__) || defined(__CYGWIN__)
#include <endian.h>
#elif defined(__APPLE__)
#include <libkern/OSByteOrder.h>
#define htobe16(x) OSSwapHostToBigInt16(x)
#define htole16(x) OSSwapHostToLittleInt16(x)
#define be16toh(x) OSSwapBigToHostInt16(x)
#define le16toh(x) OSSwapLittleToHostInt16(x)
#define htobe32(x) OSSwapHostToBigInt32(x)
#define htole32(x) OSSwapHostToLittleInt32(x)
#define be32toh(x) OSSwapBigToHostInt32(x)
#define le32toh(x) OSSwapLittleToHostInt32(x)
#define htobe64(x) OSSwapHostToBigInt64(x)
#define htole64(x) OSSwapHostToLittleInt64(x)
#define be64toh(x) OSSwapBigToHostInt64(x)
#define le64toh(x) OSSwapLittleToHostInt64(x)
#elif defined(__OpenBSD__)
#include <sys/endian.h>
#elif defined(__NetBSD__) || defined(__FreeBSD__) || defined(__DragonFly__)
#include <sys/endian.h>
#define be16toh(x) betoh16(x)
#define le16toh(x) letoh16(x)
#define be32toh(x) betoh32(x)
#define le32toh(x) letoh32(x)
#define be64toh(x) betoh64(x)
#define le64toh(x) letoh64(x)
#endif

#include <stdint.h>
#include <stdlib.h>

template <int n>
inline void swap_bytes(void* bytes) {
  char* cbytes = reinterpret_cast<char*>(bytes);
  for (int i = 0, j = n - 1; i < n / 2; ++i, --j) {
    const char tmp = cbytes[i];
    cbytes[i] = cbytes[j];
    cbytes[j] = tmp;
  }
}

template <>
inline void swap_bytes<4>(void* bytes) {
  char* cbytes = reinterpret_cast<char*>(bytes);
  const char tmp[2] = {cbytes[0], cbytes[1]};
  cbytes[0] = cbytes[3];
  cbytes[1] = cbytes[2];
  cbytes[2] = tmp[1];
  cbytes[3] = tmp[0];
}

// Determine whether the machine is bigendian or not
// NOTE: The output of this function is known at compile time, a good
// compiler should know this (i.e. GCC does know it at compile time).
inline bool is_big_endian() {
  union {
    uint32_t i;
    char c[4];
  } bint = {0x01020304};
  return bint.c[0] == 1;
}

// Convert float number in big-endian to the host endianness
// NOTE: A clever compiler should optimize this heavily
inline float beftoh(float f) {
  if (is_big_endian()) return f;
  float f2 = f;
  swap_bytes<sizeof(float)>(&f2);
  return f2;
}

// Convert float in host endianness to big-endian
// NOTE: A clever compiler should optimize this heavily
inline float htobef(float f) {
  if (is_big_endian()) return f;
  float f2 = f;
  swap_bytes<sizeof(float)>(&f2);
  return f2;
}

#endif  // FAST_PCA_ENDIAN_H_
