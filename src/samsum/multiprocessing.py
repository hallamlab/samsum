import queue
from itertools import count
from threading import Thread
from multiprocessing import Process, Queue


class MapWorker(Process):
    """
    Process that reads items from an input_queue, applies a func to them and puts them on an output_queue
    """
    def __init__(self, func, input_queue, output_queue):
        super().__init__()
        self.func = func
        self.input_queue = input_queue
        self.output_queue = output_queue

    def run(self):
        while True:
            item = self.input_queue.get()
            if item is StopIteration:
                break
            k, v = item
            self.output_queue.put((k, self.func(v)))


class ProcessMap(Thread):

    def __init__(self, func, iterator, n_proc, output_queue=None):
        super().__init__()
        self.key_map = {}
        self.iterator = iterator
        self.work_queue = Queue(n_proc * 2)
        self.output_queue = output_queue or Queue()
        self.processes = [MapWorker(func, self.work_queue, self.output_queue) for _ in range(n_proc)]

    def start(self):
        for process in self.processes:
            process.start()
        super().start()

    def run(self):
        for (k, v) in self.iterator:
            self.work_queue.put((id(k), v))
            self.key_map[id(k)] = k
        for _ in self.processes:
            self.work_queue.put(StopIteration)
        for process in self.processes:
            process.join()
        self.output_queue.put(StopIteration)

    def __iter__(self):
        self.start()
        while True:
            item = self.output_queue.get()
            if item is StopIteration:
                break
            k, v = item
            yield self.key_map.pop(k), v


class MapWorkerThread(Thread):
    """
    Process that reads items from an input_queue, applies a func to them and puts them on an output_queue
    """
    def __init__(self, func, input_queue=None, output_queue=None):
        super().__init__()
        self.func = func
        self.input_queue = input_queue
        self.output_queue = output_queue

    def run(self):
        while True:
            item = self.input_queue.get()
            if item is StopIteration:
                self.output_queue.put(item)
                break
            k, v = item
            self.output_queue.put((k, self.func(v)))


class ThreadMap(Thread):

    def __init__(self, worker_type, iterator, n_thread, maxsize=2):
        super().__init__()
        self.iterator = iterator
        self.n_thread = n_thread
        self.work_queues = [queue.Queue(maxsize) for _ in range(n_thread)]
        self.output_queues = [queue.Queue(maxsize) for _ in range(n_thread)]
        self.workers = [worker_type(input_queue=in_q, output_queue=out_q) for (in_q, out_q) in zip(self.work_queues, self.output_queues)]

    def start(self):
        for worker in self.workers:
            worker.start()
        super().start()

    def __iter__(self):
        self.start()
        for i in count():
            item = self.output_queues[i % self.n_thread].get()
            if item is StopIteration:
                #do we need to empty output_queues in order to join worker threads?
                for j in range(i + 1, i + self.n_thread):
                    self.output_queues[j % self.n_thread].get()
                break
            yield item

    def run(self):
        for i, (k, v) in enumerate(self.iterator):
            self.work_queues[i % self.n_thread].put((k, v))
        for q in self.work_queues:
            q.put(StopIteration)
        for worker in self.workers:
            worker.join()


