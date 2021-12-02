import queue
import tqdm
import multiprocessing
from itertools import count
from threading import Thread


class ThreadMap(Thread):
    ###############################################################################
    #                                                                             #
    #    This class was developed by Oxford Nanopore Technologies (ONT) as part   #
    #    of the software bonito: https://github.com/nanoporetech/bonito.          #
    #                                                                             #
    ###############################################################################
    def __init__(self, worker_type, iterator, n_thread, maxsize=2):
        super().__init__()
        self.iterator = iterator
        self.n_thread = n_thread
        self.work_queues = [queue.Queue(maxsize) for _ in range(n_thread)]
        self.output_queues = [queue.Queue(maxsize) for _ in range(n_thread)]
        self.workers = [worker_type(input_queue=in_q, output_queue=out_q) for (in_q, out_q) in
                        zip(self.work_queues, self.output_queues)]

    def start(self):
        for worker in self.workers:
            worker.start()
        super().start()

    def __iter__(self):
        self.start()
        for i in count():
            item = self.output_queues[i % self.n_thread].get()
            if item is StopIteration:
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


def tqdm_multiprocessing(func, arguments_list, num_processes: int, pbar_desc: str, disable=False) -> list:
    n_tasks = len(arguments_list)

    if n_tasks == 0:
        return []
    pool = multiprocessing.Pool(processes=num_processes)

    jobs = []
    result_list_tqdm = []
    pbar = tqdm.tqdm(jobs, total=n_tasks, desc=pbar_desc, ncols=120, disable=disable)

    def update(*a):
        pbar.update()

    for args in arguments_list:
        jobs.append(pool.apply_async(func=func, args=(*args,), callback=update))
    pool.close()

    for job in pbar:
        result_list_tqdm.append(job.get())

    pbar.close()

    return result_list_tqdm
