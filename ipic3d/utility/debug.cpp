/* This release was prepared by Dana Akhmetova <danaak@kth.se>/<danieka@gmail.com> on behalf of the INTERTWinE European Exascale Project <http://www.intertwine-project.eu> */
/* iPIC3D was originally developed by Stefano Markidis and Giovanni Lapenta. 
 * This release was contributed by Alec Johnson and Ivy Bo Peng.
 * Publications that use results from iPIC3D need to properly cite  
 * 'S. Markidis, G. Lapenta, and Rizwan-uddin. "Multi-scale simulations of 
 * plasma with iPIC3D." Mathematics and Computers in Simulation 80.7 (2010): 1509-1519.'
 *
 *        Copyright 2015 KTH Royal Institute of Technology
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at 
 *
 *         http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


#ifndef NO_MPI
  #include "MPIdata.h" // for get_rank
#endif
#include "ompdefs.h" // for omp_get_thread_num
#include "debug.h"

#define implement_dprintvar_fileLine(code,type) \
  void printvar_fileLine(const char* func, const char* file, int line, \
    const char* name, type val) \
  { \
    fprintf_fileLine(stdout,"DEBUG", func,file,line, code " == %s",val,name); \
  }

implement_dprintvar_fileLine("%s", const char *);
implement_dprintvar_fileLine("%d", int);
implement_dprintvar_fileLine("%lld", long long);
implement_dprintvar_fileLine("%e", double);
implement_dprintvar_fileLine("%p", const void *);

// void dfprintf_fileLine(FILE * fptr, const char *func, const char *file, int line_number, const char *format, ...)
// {
//   // writing directly to fptr would avoid limiting the length
//   // of the output string, but by first writing to a string
//   // we achieve thread safety.
//   //
//   // write the message to a string.
//   //
//   const int maxchars = 1024;
//   char error_msg[maxchars+2];
//   // identify the process and thread
//   char process_thread_str[20];
//   #ifndef NO_MPI
//     #ifdef _OPENMP
//       snprintf(process_thread_str, 20, "(%d.%d) ",
//         MPIdata::get_rank(), omp_get_thread_num());
//     #else
//       snprintf(process_thread_str, 20, "(%d)",
//         MPIdata::get_rank());
//     #endif
//   #else
//     #ifdef _OPENMP
//       snprintf(process_thread_str, 20, "(.%d) ",
//         omp_get_thread_num());
//     #else
//       snprintf(process_thread_str, 20, "");
//     #endif
//   #endif
//   char *sptr = error_msg;
//   int chars_so_far=0;
//   va_list args;
//   va_start(args, format);
//   chars_so_far = snprintf(sptr, maxchars,
//     "%sDEBUG %s(), %s:%d: ",
//     process_thread_str,
//     func, file, // my_basename(file),
//     line_number);
//   /* print out remainder of message */
//   chars_so_far += vsnprintf(sptr+chars_so_far, maxchars-chars_so_far, format, args);
//   va_end(args);
//   sprintf(sptr+chars_so_far, "\n");
// 
//   // print the message
//   fflush(fptr);
//     fprintf(fptr,error_msg);
//   fflush(fptr);
// }

// implemented here because the new
// standard basename() in libgen.h eliminates the const
// (because standard basename returns non-const? -- seems
// like an attempted correction in the wrong direction...).
const char *my_basename (const char *name)
{
  const char *base;
  for (base = name; *name; name++)
    if (*name == '/') base = name + 1;
  return base;
}

void fprintf_fileLine(FILE * fptr,
  const char *type, const char *func, const char *file, int line_number,
  const char *format, ...)
{
  //if(MPIdata::get_rank()) return;
  // writing directly to fptr would avoid limiting the length
  // of the output string, but by first writing to a string
  // we achieve thread safety.
  //
  // write the message to a string.
  //
  const int maxchars = 1024;
  char error_msg[maxchars+2];
  // identify the process and thread
  char process_thread_str[20];
  #ifndef NO_MPI
    #ifdef _OPENMP
      snprintf(process_thread_str, 20, "(%d.%d) ",
        MPIdata::get_rank(), omp_get_thread_num());
    #else
      snprintf(process_thread_str, 20, "(%d)",
        MPIdata::get_rank());
    #endif
  #else
    #ifdef _OPENMP
      snprintf(process_thread_str, 20, "(.%d) ",
        omp_get_thread_num());
    #else
      snprintf(process_thread_str, 20, "");
    #endif
  #endif
  char *sptr = error_msg;
  int chars_so_far=0;
  va_list args;
  va_start(args, format);
  chars_so_far = snprintf(sptr, maxchars,
    "%s%s %s(), %s:%d: ",
    process_thread_str,
    type,
    func, /*file,*/ my_basename(file),
    line_number);
  /* print out remainder of message */
  chars_so_far += vsnprintf(sptr+chars_so_far, maxchars-chars_so_far, format, args);
  va_end(args);
  sprintf(sptr+chars_so_far, "\n");

  // print the message
  fflush(fptr);
    fprintf(fptr,error_msg);
  fflush(fptr);
}

