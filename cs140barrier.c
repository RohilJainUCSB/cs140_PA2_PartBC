/*
 * File: cd140barrier.c
 *
 */

#include "cs140barrier.h"

/******************************************************
 * Initialize variables in a cs140barrier.
 * It should initize the assoicated phtread mutex and condition variables
 * It should also set the intial values for other fields.
 * The odd_rount field should be intiaized to be False.
 *
 * Argument:  bstate -- keep the state of a cs140barrier,
 *                      where bstate->odd_round is initialized to be False.
 *            total_nthread -- the number of threads using this barrier.
 *
 * Return:   0 successful, otherwise -1 meaning failed.
 */

int cs140barrier_init(cs140barrier *bstate, int total_nthread) {
    //Handle error cases where there is no barrier value or no threads or negative threads
    if(bstate == NULL || total_nthread <= 0)
    {
        return -1;
    }

    if(pthread_mutex_init(&bstate->barrier_mutex, NULL) != 0) //Initialize the mute
    {
        return -1;
    }
    if(pthread_cond_init(&bstate->barrier_cond, NULL) != 0) //Initialize the condition var
    {
        return -1;
    }

    bstate->total_nthread = total_nthread;
    bstate->arrive_nthread = 0; //No threads have arrived yet
    bstate->odd_round = False;   //Start with False (even round)

    return 0;
}

/******************************************************
 * This function shall synchronize participating threads at the barrier
 * referenced by bstate. The calling thread shall block
 * until the required number of threads have called cs140barrier_wait().
 * This completes this round of synchronization.
 *
 * Argument:  bstate -- keep the state of a cs140barrier.
 *
 * Return:   Upon successful completion, the last thread which observes all
 *           threads have arrived in this barrier round should return 1,
 *           and other threads should return 0.
 *
 * Algorithm: Use the sense reversal technique discussed in class.
 *            If a thread that calls this function is not the last thread,
 *            it should conditionally wait until the round flag changes.
 *            The last thread should change the flag and other fields.
 */

int cs140barrier_wait(cs140barrier *bstate) {
    int returnValue = 0;

    pthread_mutex_lock(&bstate->barrier_mutex);

    boolean my_round = bstate->odd_round;
    bstate->arrive_nthread++; //Increase the arrived threads

    //If this is the last thread to arrive
    if(bstate->arrive_nthread == bstate->total_nthread)
    {
        //Reset the arrived threads
        bstate->arrive_nthread = 0;
        //Flip the sense attribute and broadcast the global sense
        bstate->odd_round = (bstate->odd_round == False) ? True : False;
        pthread_cond_broadcast(&bstate->barrier_cond);
        returnValue = 1; //Return value for the last thread is 1
    }
    else
    {
        //Wait until odd_round changes 
        while(bstate->odd_round == my_round)
        {
            pthread_cond_wait(&bstate->barrier_cond, &bstate->barrier_mutex);
        }
        returnValue = 0;
    }

    pthread_mutex_unlock(&bstate->barrier_mutex);
    return returnValue;
}

/******************************************************
 * Destroy mutex and cond variables in a cs140barrier.
 * Note that the memory of bstate is not freed here.
 *
 * Argument:  bstate -- keep the state of a cs140barrier,
 *
 * Return:   0 successful, otherwise an error number is returned.
 */

int cs140barrier_destroy(cs140barrier *bstate) {
    if(pthread_mutex_destroy(&bstate->barrier_mutex) != 0)
    {
        return -1;
    }
    if(pthread_cond_destroy(&bstate->barrier_cond) != 0)
    {
        return -1;
    }

    return 0;
}
