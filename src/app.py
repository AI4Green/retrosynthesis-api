from multiprocessing import Process

from sources import app
from sources.retrosynthesis.worker import results, queue, worker


if __name__ == '__main__':
	# Set up the worker process in the background to run
	# the retrosynthesis without blocking the API
	p = Process(target=worker, args=(queue, results), daemon=True)
	p.start()

	# Run the API
	app.run(debug=False, host='0.0.0.0', port=8000)
