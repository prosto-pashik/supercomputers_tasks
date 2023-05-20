#include <cstdlib>
#include <iostream>
#include <cstring>
#include <pthread.h>
#include <ctime>

using namespace std;

const int NUM_THREADS = 5;

struct args{
    int id;
    const char *msg;
    int out;
};

void* thread_job(void* arg) {
//    pthread_attr_t* attr = (pthread_attr_t *) arg;
    args* pArguments = (args *) arg;

    int x = 0;
    clock_t local_start = clock();
    for(int i=0; i<10000000; ++i) {
        ++x;
    }
    clock_t local_end = clock();
//    cout << "Run operations time: " << (double)(local_end - local_start) / CLOCKS_PER_SEC << endl;
    std::cout << "Compute for " << pArguments->id << " thread with message " << pArguments->msg << " \n";

    pthread_exit(nullptr);
    return ((void*) x);
}


void *thread_job1(void *arg)
{
    int *id = (int*) arg; // получаем идентификатор потока
//    cout << "Thread " << *id << endl;

    int x = 0;
    clock_t local_start = clock();
    for(int i=0; i<10000000; i++) {
        ++x;
    }
    clock_t local_end = clock();
//    cout << "Run operations time: " << (double)(local_end - local_start) / CLOCKS_PER_SEC << endl;

    pthread_exit(nullptr);
}

int main() {
    pthread_t threads[NUM_THREADS];
    pthread_attr_t thread_attrs[NUM_THREADS];
    const char *messages[] = {"One", "Two", "Three", "Four", "Five"};

    args args[NUM_THREADS];
    for (int i = 0; i < NUM_THREADS; ++i) {
        args[i].id = i + 1;
        args[i].msg = messages[i];
    }

    int err;
    int ids[NUM_THREADS];

    for(int i = 0; i < NUM_THREADS; i++) {
        err = pthread_attr_init(&thread_attrs[i]);
        if (err != 0) {
            cout << "Cannot create thread attribute: " << strerror(err) << endl;
            exit(-1);
        }
        pthread_attr_setdetachstate(&thread_attrs[i], PTHREAD_CREATE_JOINABLE);
        pthread_attr_setstacksize(&thread_attrs[i], 1024 * 1024 * i);


        ids[i] = i + 1;
        clock_t start = clock();
        err = pthread_create(&threads[i], &thread_attrs[i], thread_job, (void *) &args[i]);

        if(err != 0) {
            cout << "Cannot create a thread: " << strerror(err) << endl;
            exit(-1);
        }
        clock_t end = clock();
        double time_per_thread = (double)(end - start) / CLOCKS_PER_SEC;
//        cout << "Create thread time: " << time_per_thread << endl;
    }

    for (auto & thread : threads) {
        pthread_join(thread, nullptr);
    }

    for (auto & thread_attr : thread_attrs) {
        pthread_attr_destroy(&thread_attr);
    }


    pthread_exit(nullptr);
}