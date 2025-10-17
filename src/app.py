import logging
from multiprocessing import Process

from sources import app
from sources.retrosynthesis.worker import get_shared_objects, worker


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(processName)s] %(levelname)s: %(message)s",
    )
    # Set up shared queue and results
    from sources.retrosynthesis import routes

    queue, results = get_shared_objects()
    routes.queue = queue
    routes.results = results

    # Set up the worker process in the background to run
    # the retrosynthesis without blocking the API
    logging.info("Starting retrosynthesis sub-process...")
    p = Process(target=worker, args=(queue, results), daemon=False)
    p.start()

    # Run the API
    logging.info("Starting retrosynthesis API...")
    app.run(debug=False, host="0.0.0.0", port=8000)
