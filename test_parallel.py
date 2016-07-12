from multiprocessing import Process, Lock
import os

# def f(l, i):
#     l.acquire()
#     try:
#         print('hello world', i)
#     finally:
#         l.release()



def f(l, i):
    print("power: ",l**i)
    

if __name__ == '__main__':
    lock = 2
    print("main process: ",os.getpid())

    for num in range(10):
        Process(target=f, args=(lock, num)).start()

# lock = 2
# print("main process: ",os.getpid())

# for num in range(10):
#     Process(target=f, args=(lock, num)).start()