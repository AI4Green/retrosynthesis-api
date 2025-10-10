from multiprocessing import Process, Queue, Manager
import time, uuid, os

_manager = Manager()
results = _manager.dict()
_queue = Queue()