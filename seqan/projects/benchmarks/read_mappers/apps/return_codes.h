/*==========================================================================
   SeqAn - The Library for Sequence Analysis
   http://www.seqan.de 
  ==========================================================================
   Copyright (C) 2010
  
   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 3 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   Lesser General Public License for more details.
  
  ==========================================================================
   Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
  ==========================================================================
   Return codes for the read mapper benchmark tools.
  ==========================================================================*/

#ifndef BENCHMARKS_READ_MAPPERS_RETURN_CODES_H_
#define BENCHMARKS_READ_MAPPERS_RETURN_CODES_H_

// Define some return codes.
const int kRetOk = 0;       // OK, no errors.
const int kRetArgsErr = 1;  // Errors in arguments.
const int kRetIoErr = 2;    // I/O error, problem reading files.
const int kFatalErr = 3;    // Some other sort of fatal error.

#endif  // BENCHMARKS_READ_MAPPERS_RETURN_CODES_H_
