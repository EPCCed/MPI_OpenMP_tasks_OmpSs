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

#ifndef __TimeTasks_H__
#define __TimeTasks_H__
#include "assert.h"
#include "errors.h"

/* Avoid direct use of this class.
   Instead, use and add to the macros at the bottom
   so that we can redefine the macros when desired
   (e.g. defining them to the empty string to
   remove performance penalty).
 */

class TimeTasks
{
 public:

  // legitimate active subcycle values
  //
  // timeTasks_set_task(0) is a no-op, so
  // MOMENT_REDUCTION=0
  // would prevent monitoring of this task.
  //
  enum Tasks // order must agree with taskNames in TimeTasks.cpp
  {
    NONE = 0,
    //
    FIELDS,
    PARTICLES,
    MOMENTS,
    after_exclusive,
    //
    COMMUNICATING,
    before_communication,
    FLDS_COMM,
    FLDS_MPI_ALLREDUCE,
    FLDS_MPI_SENDRECV,
    PCLS_COMM,
    PCLS_MPI_ALLREDUCE,
    PCLS_MPI_SENDRECV,
    MOMS_COMM,
    MOMS_MPI_ALLREDUCE,
    MOMS_MPI_SENDRECV,
    //
    before_report_list,
    REDUCE_FIELDS,
    BFIELD,
    MOMENT_PCL_SORTING,
    MOMENT_ACCUMULATION,
    MOMENT_REDUCTION,
    MOVER_PCL_SORTING,
    MOVER_PCL_MOVING,
    TRANSPOSE_PCLS_TO_AOS,
    TRANSPOSE_PCLS_TO_SOA,
    WRITE_FIELDS,
    WRITE_PARTICLES,
    PCLS_MPI_Isend,
    PCLS_MPI_Irecv,
    PCLS_MPI_Wait,
    PCLS_MPI_Cancel,
    PCLS_MPI_Request_free,
    PCLS_MPI_Test,
    PCLS_MPI_Waitany,
    //
    NUMBER_OF_TASKS // this line should be last
  };

 private:
  //enum Modes // for exclusive tasks
  //{
  //  COMPUTATION = 0,
  //  COMMUNICATION,
  //};

 public: // methods

  TimeTasks() {
    resetCycle();
  }

  // monitoring
  //
  void resetCycle();
  //
  // hack to support averaging timeTasks copies of all threads.
  //
  void operator+=(const TimeTasks& arg);
  void operator/=(int num);
  void operator=(const TimeTasks& arg);
  //
  // provide start_time on ending call
  //
  //void end_communicating(double start_time);
  void end_sendrecv(double start_time);
  void end_allreduce(double start_time);
  void start_main_task(TimeTasks::Tasks taskid);
  void end_main_task(TimeTasks::Tasks taskid, double start_time);
  void start_task(TimeTasks::Tasks taskid);
  void end_task(TimeTasks::Tasks taskid, double start_time);
  //
  // provide start_time at starting call
  //
  void start_task(TimeTasks::Tasks taskid, double start_time);
  void end_task(TimeTasks::Tasks taskid);

  // accessors
  //
  bool is_active(Tasks taskid){
    bool retval = active[taskid];
    //if(retval&& stack_depth[taskid]==0)
    //{
    //  eprintf("active task %s has depth %d",
    //    get_taskname(taskid), stack_depth[taskid]);
    //}
    return retval;
  }
  //bool get_communicating() { return communicating; }
  //void set_communicating(bool val) { communicating = val; }
  int get_stack_depth(TimeTasks::Tasks taskid) { return stack_depth[taskid]; }

  // reporting
  //
 private:
  void print_cycle_times(int cycle, double* task_duration,
    const char* reduce_mode="avg");
  void print_cycle_times(int cycle, const char* reduce_mode);
 public:
  void print_cycle_times(int cycle);

 private:

  // is task exclusive?
  bool is_exclusive(Tasks taskid) { return (taskid < after_exclusive); }

  // reporting
  //
  //double get_time(int arg) {
  //  return task_duration[arg];
  //}
  //double get_communicate(int arg) {
  //  return communicate[arg];
  //}
  //double get_compute(int arg) {
  //  return get_time(arg) - get_communicate(arg);
  //}
  const char* get_taskname(int arg);

 private:
  int active_task;
  bool active[NUMBER_OF_TASKS];
  //bool communicating;
  double task_duration[NUMBER_OF_TASKS];
  //double communicate[NUMBER_OF_TASKS];
  //double sendrecv[NUMBER_OF_TASKS];
  //double allreduce[NUMBER_OF_TASKS];
  int stack_depth[NUMBER_OF_TASKS];
  double start_times[NUMBER_OF_TASKS];
};

extern TimeTasks timeTasks;

// construct an anonymous instance of TimeTasksCaller
class TimeTasks_caller_to_set_main_task_for_scope
{
  double start_time;
  TimeTasks::Tasks task;
 public:
  TimeTasks_caller_to_set_main_task_for_scope(TimeTasks::Tasks _task);
  ~TimeTasks_caller_to_set_main_task_for_scope();
};

class TimeTasks_caller_to_set_task_for_scope
{
  bool already_active;
  double start_time;
  TimeTasks::Tasks task;
 public:
  TimeTasks_caller_to_set_task_for_scope(TimeTasks::Tasks task_);
  ~TimeTasks_caller_to_set_task_for_scope();
};

class TimeTasks_caller_to_set_communication_mode_for_scope
{
 private:
  bool already_communicating;
  double start_time;
 public:
  TimeTasks_caller_to_set_communication_mode_for_scope();
  ~TimeTasks_caller_to_set_communication_mode_for_scope();
};

// These macros could be changed to provide file and line number
//
// We need to create nonanonymous instances so that the destructor
// will not be called until the end of the scope, so we use the preprocessor
// to generate unique names of nonanonymous instances.
//
#define timeTasks_set_main_task(task) \
  TimeTasks_caller_to_set_main_task_for_scope myFunnyInstance(task);
// unfortunately this just pastes __func__ and __LINE__ literally
//#define timeTasks_set_task(task) \
//TimeTasks_caller_to_set_task_for_scope myFunnyName##__func__##__LINE__(task);
#define timeTasks_set_task(task) \
  TimeTasks_caller_to_set_task_for_scope myFunnyName(task);
#define timeTasks_set_communicating() timeTasks_set_task(TimeTasks::COMMUNICATING);
//#define timeTasks_set_communicating() \
//  TimeTasks_caller_to_set_communication_mode_for_scope myFunnyCommunicationInstance;
//
// The scoping trick does not work if the timeTasks call needs to be conditional,
// so we also provide the ability to explicitly begin and end.
#define timeTasks_begin_task(task) if(task) timeTasks.start_task(task, MPI_Wtime());
#define timeTasks_end_task(task) if(task) timeTasks.end_task(task);
//

#endif
