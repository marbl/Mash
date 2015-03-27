#ifndef ThreadPool_h
#define ThreadPool_h

#include <pthread.h>
#include <queue>

template <class TypeInput, class TypeOutput>
class ThreadPool
{
public:
    
    ThreadPool(TypeOutput * (* functionNew)(TypeInput *), unsigned int threadCountNew);
    ~ThreadPool();
    
    bool outputAvailable() const;
    TypeOutput * popOutputWhenAvailable(); // output must be deleted by calling function
    bool running() const;
    void runWhenThreadAvailable(TypeInput * input); // thread deletes input when finished
    
private:
    
    struct OutputQueueNode
    {
        // used to preserve input order when outputting
        
        OutputQueueNode * prev;
        OutputQueueNode * next;
        
        TypeOutput * output;
        bool ready;
    };
    
    unsigned int threadCount;
    
    pthread_t * threads;
    
    static void * thread(void *);
    
    TypeOutput * (* function)(TypeInput *);
    TypeInput * inputCurrent;
    OutputQueueNode * outputQueueNodeCurrent;
    
    pthread_mutex_t * mutexInput;
    pthread_mutex_t * mutexOutput;
    
    pthread_cond_t * condInput;
    pthread_cond_t * condOutput;
    
    OutputQueueNode * outputQueueHead;
    OutputQueueNode * outputQueueTail;
    
    bool finished;
    friend void * thread(void *);
};


#include "ThreadPool.hxx"

#endif
