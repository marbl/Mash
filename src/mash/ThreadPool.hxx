#include "ThreadPool.h"
#include <stdlib.h>
#include <stdio.h>

template <class TypeInput, class TypeOutput>
ThreadPool<TypeInput, TypeOutput>::ThreadPool(TypeOutput * (* functionNew)(TypeInput *), unsigned int threadCountNew)
    :
    threadCount(threadCountNew),
    function(functionNew)
{
    mutexInput = new pthread_mutex_t();
    mutexOutput = new pthread_mutex_t();
    
    condInput = new pthread_cond_t();
    condOutput = new pthread_cond_t();
    
    pthread_mutex_init(mutexInput, NULL);
    pthread_mutex_init(mutexOutput, NULL);
    
    pthread_cond_init(condInput, NULL);
    pthread_cond_init(condOutput, NULL);
    
    inputCurrent = 0;
    
    outputQueueHead = 0;
    outputQueueTail = 0;
    
    finished = false;
    
    threads = new pthread_t[threadCount];
    
    for ( int i = 0; i < threadCount; i++ )
    {
        pthread_create(&threads[i], NULL, &ThreadPool::thread, this);
    }
}

template <class TypeInput, class TypeOutput>
ThreadPool<TypeInput, TypeOutput>::~ThreadPool()
{
    pthread_mutex_lock(mutexInput);
    finished = true;
    pthread_cond_broadcast(condInput);
    pthread_mutex_unlock(mutexInput);
    
    for ( int i = 0; i < threadCount; i++ )
    {
        pthread_join(threads[i], NULL);
    }
    
    delete [] threads;
    
    while ( outputQueueHead != 0 )
    {
        OutputQueueNode * next = outputQueueHead->next;
        delete outputQueueHead;
        outputQueueHead = next;
    }
    
    delete mutexInput;
    delete mutexOutput;
    
    delete condInput;
    delete condOutput;
}

template <class TypeInput, class TypeOutput>
bool ThreadPool<TypeInput, TypeOutput>::outputAvailable() const
{
    bool available;
    
    //std::cout << "OA locking\n";
    pthread_mutex_lock(mutexOutput);
    available = outputQueueHead != 0 && outputQueueHead->ready;
    //std::cout << "OA unlocking\n";
    pthread_mutex_unlock(mutexOutput);
    
    return available;
}

template <class TypeInput, class TypeOutput>
TypeOutput * ThreadPool<TypeInput, TypeOutput>::popOutputWhenAvailable()
{
    //std::cout << "pop locking output\n";
    pthread_mutex_lock(mutexOutput);
    
    if ( outputQueueHead == 0 )
    {
        // TODO: error?
        //std::cout << "ERROR: waiting for output when no output queued\n";
        pthread_mutex_unlock(mutexOutput);
        return 0;
    }
    
    //std::cout << "Waiting for output\n";
    while ( ! outputQueueHead->ready )
    {
        pthread_cond_wait(condOutput, mutexOutput);
    }
    
    TypeOutput * output = outputQueueHead->output;
    
    OutputQueueNode * next = outputQueueHead->next;
    
    if ( outputQueueTail == outputQueueHead )
    {
        outputQueueTail = 0;
    }
    
    delete outputQueueHead;
    outputQueueHead = next;
    //std::cout << "pop new head: " << outputQueueHead << std::endl;
    pthread_mutex_unlock(mutexOutput);
    
    return output;
}

template <class TypeInput, class TypeOutput>
void ThreadPool<TypeInput, TypeOutput>::runWhenThreadAvailable(TypeInput * input)
{
    //std::cout << "waiting for thread\n";
    pthread_mutex_lock(mutexInput);
    
    while ( inputCurrent != 0 )
    {
        pthread_cond_wait(condInput, mutexInput);
    }
    
    //std::cout << "setting input\n";
    
    inputCurrent = input;
    
    pthread_mutex_unlock(mutexInput);
    pthread_cond_broadcast(condInput);
}

template <class TypeInput, class TypeOutput>
bool ThreadPool<TypeInput, TypeOutput>::running() const
{
    bool running;
    
    pthread_mutex_lock(mutexOutput);
    running = outputQueueHead != 0;
    pthread_mutex_unlock(mutexOutput);
    
    return running;
}

template <class TypeInput, class TypeOutput>
void * ThreadPool<TypeInput, TypeOutput>::thread(void * arg)
{
    //std::cout << "Thread created" << std::endl;
    ThreadPool * threadPool = (ThreadPool *)arg;
    TypeInput * input;
    
    while ( ! threadPool->finished )
    {
        // wait for input
        //
        //std::cout << "Waiting for input\n";
        pthread_mutex_lock(threadPool->mutexInput);
        //
        while ( ! threadPool->finished && threadPool->inputCurrent == 0 )
        {
            pthread_cond_wait(threadPool->condInput, threadPool->mutexInput);
        }
        
        if ( threadPool->finished )
        {
            pthread_mutex_unlock(threadPool->mutexInput);
            return 0;
        }
        //
        //std::cout << "Taking input " << threadPool->inputCurrent << std::endl;
        input = threadPool->inputCurrent;
        threadPool->inputCurrent = 0;
        
        //std::cout << "Enqueuing output\n";
        // enqueue output while input locked (to preserve order)
        //
        OutputQueueNode * outputQueueNode = new OutputQueueNode();
        outputQueueNode->next = 0;
        outputQueueNode->ready = false;
        //
        //std::cout << "locking output\n";
        pthread_mutex_lock(threadPool->mutexOutput);
        //
        //std::cout << "output locked\n";
        if ( threadPool->outputQueueHead == 0 )
        {
            //std::cout << "output queue empty; setting head\n";
            threadPool->outputQueueHead = outputQueueNode;
        }
        //
        outputQueueNode->prev = threadPool->outputQueueTail;
        //
        if ( threadPool->outputQueueTail != 0 )
        {
            //std::cout << "setting tail\n";
            threadPool->outputQueueTail->next = outputQueueNode;
        }
        //
        threadPool->outputQueueTail = outputQueueNode;
        //
        pthread_mutex_unlock(threadPool->mutexOutput);
        pthread_mutex_unlock(threadPool->mutexInput);
        
        pthread_cond_broadcast(threadPool->condInput);
        
        //std::cout << "Running\n";
        
        // run function
        //
        outputQueueNode->output = threadPool->function(input);
        
        delete input;
        
        // signal output
        //
        //std::cout << "outputting\n";
        outputQueueNode->ready = true;
        //
        pthread_mutex_lock(threadPool->mutexOutput);
        pthread_cond_broadcast(threadPool->condOutput);
        pthread_mutex_unlock(threadPool->mutexOutput);
    }
    
    return NULL;
}
