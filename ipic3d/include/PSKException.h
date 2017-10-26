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

// File: PSKException.h
/*! \file Base exception class */

#ifndef _PSK_EXCEPTION_H_
#define _PSK_EXCEPTION_H_

#include <string>
#include <vector>
#include <stdio.h>
#include <errno.h>

#include <cstring>

namespace PSK {

  /*! */
  class Exception {
    protected:
      bool _can_recover;          // /< Flag indicating recovery possible
      std::string _err_str;       // /< Exception error message string
      std::string _fn_str;        // /< Function name throwing exception
      std::string _type_str;      // /< Exception type string, including inheritance
      std::vector < std::string > _err_msgs;  // /< stack of associated messages
      int _sys_errno;             // /< system error number from bad system call

      std::string _sys_err_str; // /< system error message from bad system call
    public:
      Exception(void):_type_str("PSK::Exception") {;
      } Exception(const std::string & err_str, const std::string fn_str = "", int sys_errno = 0)
      :  _can_recover(true), _err_str(err_str), _fn_str(fn_str), _type_str("PSK::Exception"), _sys_errno(sys_errno) {
        if (sys_errno != 0)
          _sys_err_str = strerror(sys_errno);
        errno = 0;
      } Exception(const char *c_err_str, const char *c_fn_str = "", int sys_errno = 0)
      :  _can_recover(true), _err_str(c_err_str), _fn_str(c_fn_str), _type_str("PSK::Exception"), _sys_errno(sys_errno) {
        if (sys_errno != 0)
          _sys_err_str = strerror(sys_errno);
        errno = 0;
      }

      Exception(const Exception & e):_can_recover(e._can_recover), _err_str(e._err_str), _fn_str(e._fn_str), _type_str(e._type_str), _err_msgs(e._err_msgs), _sys_errno(e._sys_errno), _sys_err_str(e._sys_err_str) {;
      }

      void append_err_str(const std::string & str) {
        _err_str += str;
      }

      void prepend_err_str(const std::string & str) {
        _err_str = str + _err_str;
      }

      void push(const std::string & err_msg) {
        _err_msgs.push_back(err_msg);
      }

      void diag_cout() {
        printf("%s %s in %s\n",
          _type_str.c_str(), _err_str.c_str(), _fn_str.c_str());
        //std::cout << _type_str << " " << _err_str << " in " << _fn_str << std::endl;
        if (_sys_errno != 0)
        {
          printf("sys errno: %d [%s]\n", _sys_errno, _sys_err_str.c_str());
          //std::cout << "sys errno: " << _sys_errno << " [" << _sys_err_str << "]" << std::endl;
        }
        int n_msgs = _err_msgs.size();
        for (int i = 0; i < n_msgs; ++i)
        {
          printf(".. %s\n", _err_msgs[i].c_str());
          //std::cout << ".. " << _err_msgs[i] << std::endl;
        }
      }

  };

}                               // end namespace PSK

#endif // _PSK_EXCEPTION_H_
