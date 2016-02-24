// Implement your solutions in this file
#include "findmotifs.h"
#include "hamming.h"
#include <mpi.h>
#include <queue>
#include <iostream>
#include <ctime>
#include <iomanip>

#define WORKTAG 1
#define TERMINTAG 2

//#ifndef MPI_WTIME_IS_GLOBAL
#define MPI_WTIME_IS_GLOBAL true 
//#endif

std::vector<bits_t> findmotifs_worker(const unsigned int n,
                       const unsigned int l,
                       const unsigned int d,
                       const bits_t* input,
                       const unsigned int startbitpos,
                       bits_t start_value)
{
    std::vector<bits_t> results;
    bits_t pre_scheme = 0;
    int bound_of_diff = d - hamming(start_value, input[0]);
    int diff = 0;

    while( diff >= 0)
    {
        bool is_sol = true;
        bits_t potential_sol = start_value ^ (pre_scheme << startbitpos);

        for(unsigned int i = 1; i < n; i++)
        {
            if(hamming(potential_sol, input[i]) > d)
            {
                is_sol = false;
                break;
            }
        }
        if (is_sol)
            results.push_back(potential_sol);
        diff = nextModiScheme(pre_scheme, diff, l-startbitpos, bound_of_diff);
    }
    return results;
}

void worker_main()
{
    MPI_Status status;

    std::deque<bits_t> subproblem;
    std::deque<std::vector<bits_t> > results;
    std::deque<MPI_Request> Isend_request;
    int rank,flag;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);


    // TODO:
    // 1.) receive input from master (including n, l, d, input, master-depth)

    unsigned int* configure = new unsigned int[4];
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(configure, 4, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    unsigned int n, d, l, master_depth;
    n = configure[0]; l = configure[1]; d = configure[2]; master_depth = configure[3];

    bits_t *input = new bits_t[n];
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(input, n, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    // 2.) while the master is sending work:
    //      a) receive subproblems 

    bits_t receive_subprob;
    while(1){
        //For every while loop
        //1.Finish one task
        //2.Send back the result
        //3.Receive possible new task
        //4.Receive possible terminate information(Only when there's no subproblem in buff

        if (subproblem.empty())
        {
            MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            //check whether this is the terminate message, break if it is the terminate tag
            if (status.MPI_TAG == TERMINTAG)
            {
                //std::cout << "The end of the slave of Processor " << rank << std::endl;
                return;
            }        
            //Else, it should be a task message
            MPI_Recv(&receive_subprob, 1, MPI_UNSIGNED_LONG, 0, WORKTAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            subproblem.push_back(receive_subprob);

            //Recheck, in case there're two tasks
            MPI_Iprobe(0, WORKTAG, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
            if(flag)
            {
                MPI_Recv(&receive_subprob, 1, MPI_UNSIGNED_LONG, 0, WORKTAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                subproblem.push_back(receive_subprob);
            }
        }
        else
        {
            //The subproblem should contain exactly one problem(Unsolved)
            //Check whether there's new assignment. 
            MPI_Iprobe(0, WORKTAG, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
            if(flag)
            {
                MPI_Recv(&receive_subprob, 1, MPI_UNSIGNED_LONG, 0, WORKTAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                subproblem.push_back(receive_subprob);
            }

        }
        if (subproblem.empty())
          std::cout << "Wrong in Worker _ Empty\n";

        //b) locally solve (using findmotifs_worker(...))

        receive_subprob=subproblem.front();
        subproblem.pop_front();

        results.push_front(std::vector<bits_t> ());
        results[0]=findmotifs_worker(n,l,d,input,master_depth,receive_subprob);

        //c) send results to master 
        Isend_request.push_front(MPI_Request ());

        int count = results[0].size();
        if (count == 0)
        {
            results[0].push_back(~input[0]);
            count++;
        }

        MPI_Isend(&results[0][0], count, MPI_UNSIGNED_LONG, 0, WORKTAG, MPI_COMM_WORLD, &Isend_request[0]);

        //d)Clear
        int flag = 0;
        MPI_Test(&Isend_request[Isend_request.size() - 1], &flag, MPI_STATUS_IGNORE);
        if (flag)
        {
            results[results.size() - 1].clear();
            results.pop_back();
            Isend_request.pop_back();
        }
    }
}



std::vector<bits_t> findmotifs_master(const unsigned int n,
                                      const unsigned int l,
                                      const unsigned int d,
                                      const bits_t* input,
                                      const unsigned int till_depth)
{
    std::vector<bits_t> results;
    bits_t master_scheme = 0;
    int diff = 0;

    while( diff >= 0)
    {
        results.push_back(input[0]^master_scheme);
        diff = nextModiScheme(master_scheme, diff, till_depth, d);
    }
    return results;
}


struct Status_Cmp{
    bool operator()(const double* a, const double*b)
    {
        //a, b is a triplet with [#task, time, rank] 
        //return the comparison result whether a < b
        if(b[0] < a[0])
          return true;
        else if (b[0] > a[0])
          return false;
        else if (b[1] < a[1])
          return true;
        else 
          return false;
    }
};

std::vector<bits_t> master_main(unsigned int n, unsigned int l, unsigned int d,
            const bits_t* input, unsigned int master_depth)
{
    // 1.) send input to all workers (including n, l, d, input, depth)

    //BroadCast basic message, n, l, d, input, depth
    unsigned int configure[] = {n, l, d, master_depth};
    bits_t* input_broadcast = new bits_t [n];
    for(unsigned int i = 0; i < n; i++)
      input_broadcast[i] = input[i];

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(configure, 4, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(input_broadcast, n, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    delete input_broadcast;

    // 2.) solve problem till depth `master_depth` and then send subproblems
    //     to the workers and receive solutions in each communication
    //     Use your implementation of `findmotifs_master(...)` here.

    std::vector<bits_t> results; //Final result
    int nproc;  MPI_Comm_size(MPI_COMM_WORLD, &nproc); //Number of processors

    //init stores all the triplets of the worker processors' working status
    std::vector<double*> init(nproc-1);

    //update_time_pro stores the most recent update time of each processor(To prevent misuse of outdated working status)
    std::vector<double> update_time_pro(nproc);

    //Initilize init and update_time_pro
    double init_time = MPI_Wtime();
    for(int i = 1; i < nproc; i++)
    {
        init[i-1] = new double[3];
        init[i-1][0] = 0;
        init[i-1][1] = init_time; 
        init[i-1][2] = i;
        update_time_pro[i] = init_time;
    }

    //assignment_manager is a priority queue storing the status of each processor[#task, update_time, rank],
    //the idlest one will be placed at the top(By the reloaded () operator(Status_Cmp))
    //choose the idlest processor to assign new task;
    std::priority_queue<double*, std::vector<double*>, Status_Cmp> assignment_manager(init.begin(), init.end());
    std::deque<bits_t> master_solution;
    std::deque<MPI_Request> Isend_request;
    std::vector<int> taskNum(nproc);  taskNum[0] = 0;


    MPI_Status status;
    MPI_Request request; 
    int flag;
    double* proc_status;
    bits_t Scheme = 0;
    int diff = 0;

    //Start Assigning Tasks master_results[i]
    while ( diff >= 0)
    {
        if (assignment_manager.empty())
          std::cout << "Something wrong with the assignment manager\n";

        //Choose the idlest processor
        proc_status = assignment_manager.top();
        assignment_manager.pop();

        //First check whether this information is up-to-date information
        if(proc_status[1] >= update_time_pro[(int)proc_status[2]])
        {
            //New Task Partial Solution
            master_solution.push_front(Scheme ^ input[0]);
            Isend_request.push_front(request);

            //Update  scheme
            diff = nextModiScheme(Scheme, diff, master_depth, d);

  //        printf("rank %d time %.3f task %d presenttime %.3f\n", (int)proc_status[2], proc_status[1], (int)proc_status[0], MPI_Wtime());
            //Updated Time && Assign new task;
            MPI_Isend(&master_solution[0], 1, MPI_UNSIGNED_LONG, (int)proc_status[2], WORKTAG, MPI_COMM_WORLD, &Isend_request[0]);
            update_time_pro[(int)proc_status[2]] = MPI_Wtime();
            taskNum[(int)proc_status[2]] += 1;

            //Check the earliest Isend Request and refresh the master_solution & Isend_request
            MPI_Test(&Isend_request[Isend_request.size()-1], &flag, MPI_STATUS_IGNORE);
            if(flag)
            {
                master_solution.pop_back();
                Isend_request.pop_back();
            }

            //Number of assignement in the asignment_manger should be smaller than 2 ==> Assign a new task for the processor
            if (taskNum[(int)proc_status[2]] < 2)
            {
                proc_status[0] = taskNum[(int)proc_status[2]];
                proc_status[1] = update_time_pro[(int)proc_status[2]];
                assignment_manager.push(proc_status);
            }
            else
              //After assignment, this processor has 2 tasks in local memory
              //No need to put the updated status back to PQ
              delete proc_status;
        }
        else
          //Out-of-Date status data
          delete proc_status;

        //Check whether all processors have 2 task by assignment_manager.empty();
        if (assignment_manager.empty())
          //Wait new message come from any workers and update queue
          MPI_Probe(MPI_ANY_SOURCE, WORKTAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        //Receive the new message including
        //1. Solution from the processor: WORKTAG
        //2. Gather all the solutions.
        while (1)
        {
            MPI_Iprobe(MPI_ANY_SOURCE, WORKTAG, MPI_COMM_WORLD, &flag, &status);
            if (flag == 0)
              break;
            //Receive solution.
            int count;
            MPI_Get_count(&status, MPI_UNSIGNED_LONG, &count);
            bits_t* proc_sol = new bits_t[count];
            MPI_Recv(proc_sol, count, MPI_UNSIGNED_LONG, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            //Update Status Information
            taskNum[status.MPI_SOURCE] -= 1 ;
            double* tmp_status = new double[3];
            tmp_status[0] = taskNum[status.MPI_SOURCE];
            tmp_status[1] = MPI_Wtime();
            tmp_status[2] = status.MPI_SOURCE;
            assignment_manager.push(tmp_status);
            update_time_pro[status.MPI_SOURCE] = tmp_status[1];

            if (proc_sol[0] == ~input[0])
              //When no solution found, the processor will return ~input[0]
              continue;

            //Merge the results.
            results.insert(results.end(), proc_sol, proc_sol + count);
            delete proc_sol;
        }
    }

    // 3.) receive last round of solutions
    int TERMINATE = 0;

    //First, Empty the Assignment Manager's memory
    while (!assignment_manager.empty())
    {
        double* tmp_status = assignment_manager.top();
        assignment_manager.pop();
        delete tmp_status;
    }

    int task_remain = 0 ;
    for(int k = 1; k < nproc; k++)
    {
        if (taskNum[k] == 0)
          MPI_Send(&TERMINATE, 1, MPI_INT, k, TERMINTAG, MPI_COMM_WORLD);
        else
          task_remain += taskNum[k];
    }

    //Second, Collect the last round of solutions
    int task_count = 0;
    while (task_count < task_remain)
    {
        //Wait for new information(The information is sent when one subwork has been finished in workers
        MPI_Probe(MPI_ANY_SOURCE, WORKTAG, MPI_COMM_WORLD, &status);

        int count;
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &count);
        bits_t* proc_sol = new bits_t[count];
        MPI_Recv(proc_sol, count, MPI_UNSIGNED_LONG, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        taskNum[status.MPI_SOURCE] -= 1;

        if (taskNum[status.MPI_SOURCE] == 0)
            MPI_Send(&TERMINATE, 1, MPI_INT, status.MPI_SOURCE , TERMINTAG, MPI_COMM_WORLD);

        //Count finished job
        task_count++;

        //When no solution found, the processor will return ~input[0]
        if (proc_sol[0] == ~input[0]) continue;

        //Merge the results.
        results.insert(results.end(), proc_sol, proc_sol + count);
        delete proc_sol;
    }


    // 4.) terminate (and let the workers know)
    return results;
}

