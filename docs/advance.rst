Multiple processes
==================

Pyfastx can be used with `multiprocessing <https://docs.python.org/3.7/library/multiprocessing.html>`_ module to speed up the random access. Prior to reading sequences from subprocesses, you have to ensure that index file has been created from main process. The index file of pyfastx is just a SQLite3 database file which supports for concurrency. We provide two simple examples for using pyfastx with multiprocessing pool.

Example one
-----------

.. code:: python

    import random
    import pyfastx
    import multiprocessing as mp

    # process worker
    # randomly fetch five sequences and print to stdout
    def worker(woker_num, seq_counts):
        #recreate the Fasta object in subprocess
        fa = pyfastx.Fasta('test.fa')

        for i in random.sample(range(seq_counts), 5):
            print("worker {} print:\n{}".format(worker_num, fa[i].raw))

    if __name__ == '__main__':
        #ensure index file has been created in main process
        fa = pyfastx.Fasta('test.fa')

        #get sequence counts
        c = len(fa)

        #start the process pool
        pool = mp.Pool()

        #add five task workers to run
        for n in range(5):
            pool.apply_async(worker, args=(n, c))

        #wait for tasks to finish
        pool.close()
        pool.join()
