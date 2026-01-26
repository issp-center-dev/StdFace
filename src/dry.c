/**
 * @file dry.c
 * @brief Main entry point for the StdFace standard-input-file generator
 *
 * @details
 * This file provides the command-line driver that parses arguments and
 * invokes StdFace_main() to process an input file.  It also supports
 * a @c -v flag to print the version string.
 *
 * @copyright
 * mVMC - A numerical solver package for a wide range of quantum lattice
 * models based on many-variable Variational Monte Carlo method
 * Copyright (C) 2016 The University of Tokyo, All rights reserved.
 *
 * This program is developed based on the mVMC-mini program
 * (https://github.com/fiber-miniapp/mVMC-mini)
 * which follows "The BSD 3-Clause License".
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see http://www.gnu.org/licenses/.
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "version.h"

void StdFace_main(char *fname);

/**
 * @brief Print a usage message showing the expected command-line syntax
 *
 * @param[in] prog  Program name (typically argv[0])
 */
void usage(const char *prog)
{
  printf("usage: %s stan.in\n", prog);
}

/**
 * @brief Entry point for the StdFace command-line driver
 *
 * @details
 * Parses command-line arguments and dispatches to the appropriate action:
 * - No arguments: prints version and usage information, then exits with
 *   a non-zero status.
 * - @c -v flag: prints the version string.
 * - Otherwise: treats the first argument as an input filename and passes
 *   it to StdFace_main() for processing.
 *
 * @param[in] argc  Number of command-line arguments
 * @param[in] argv  Array of command-line argument strings
 * @return 0 on success, 1 if no input file was provided
 */
int main(int argc, char *argv[])
{
  /* Return status code indicating success/failure */
  int status = 0;

  if (argc < 2) {
    /* Print version and usage if no arguments provided */
    printVersion();
    usage(argv[0]); 
    status = 1;
  } else if (strcmp(argv[1], "-v") == 0) {
    /* Print version if -v flag specified */
    printVersion();
  } else {
    /* Process input file */
    StdFace_main(argv[1]);
  }

  return status;
}
