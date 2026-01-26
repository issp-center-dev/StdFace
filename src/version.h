/**
 * @file version.h
 * @brief Version information and printing utility for StdFace
 *
 * @copyright
 * HPhi-mVMC-StdFace - Common input generator
 * Copyright (C) 2015 The University of Tokyo
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */
#ifndef _INCLUDE_VERSION
#define _INCLUDE_VERSION

#include <stdio.h>
#include <string.h>

/* Semantic Versioning http://semver.org */
/* <major>.<minor>.<patch>-<prerelease> */
#define VERSION_MAJOR 0
#define VERSION_MINOR 5
#define VERSION_PATCH 0
#define VERSION_PRERELEASE "" /* "alpha", "beta.1", etc. */

/**
 * @brief Print the StdFace version string to standard output
 */
void printVersion() {
  printf("StdFace version %d.%d.%d", VERSION_MAJOR, VERSION_MINOR,
         VERSION_PATCH);
  if (strlen(VERSION_PRERELEASE) > 0) {
    printf("-%s", VERSION_PRERELEASE);
  }
  printf("\n");
  return;
};

#endif /* _INCLUDE_VERSIN */
