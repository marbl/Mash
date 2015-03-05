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
    
    unsigned int threadCount;
    
    pthread_t * threads;
    
    static void * thread(void *);
    
public: // TODO: these need to be accessed by static thread(). Better way?
    
    struct OutputQueueNode
    {
        // used to preserve input order when outputting
        
        OutputQueueNode * prev;
        OutputQueueNode * next;
        
        TypeOutput * output;
        bool ready;
    };
    
    TypeOutput * (* function)(TypeInput *);
    TypeInput * inputCurrent;
    
    pthread_mutex_t * mutexInput;
    pthread_mutex_t * mutexOutput;
    
    pthread_cond_t * condInput;
    pthread_cond_t * condOutput;
    
    OutputQueueNode * outputQueueHead;
    OutputQueueNode * outputQueueTail;
    
    bool finished;
};


#include "ThreadPool.hxx"

#endif
