import os
import json
import time
import uuid
from flask import request, jsonify, current_app

import sources.retrosynthesis.startup
from sources import app
from sources.retrosynthesis.worker import queue, results

# ACCESS_KEYS can be a comma-separated string of keys per client
# e.g. client1-key,client2-key
ACCESS_KEYS = os.getenv("KEYS", "retro_key").split(",")


@app.route("/", methods=["GET"])
def service_check():
    finder = sources.retrosynthesis.startup.make_config()
    # stock1 = finder.config.stock
    page_data = {
        "Message": "Retrosynthesis service is running",
        "Timestamp": time.time(),
        "Version": "v2.05.25",
    }
    json_dump = json.dumps(page_data)
    return json_dump


@app.route("/retrosynthesis_api/", methods=["GET"])
def retrosynthesis():
    access_key = str(request.args.get("key"))
    if access_key not in ACCESS_KEYS:
        print("Invalid key")
        return json.dumps({"Message": "Invalid key", "Timestamp": time.time()})

    # Get job parameters
    smiles = str(request.args.get("smiles"))
    enhancement = str(request.args.get("enhancement", "Default"))
    iteration_limit = int(request.args.get("iterations"))
    max_transforms = int(request.args.get("max_depth"))
    time_limit = int(request.args.get("time_limit"))

    # Create unique job ID
    job_id = str(uuid.uuid4())

    # Queue the analysis
    current_app.logger.info(f"Queueing job {job_id}")
    queue.put(
        (job_id, smiles, enhancement, iteration_limit, max_transforms, time_limit)
    )

    return jsonify({"job_id": job_id}), 200


@app.route("/results/<job_id>")
def get_results(job_id: str):
    """Get the results for a retrosynthesis job with the given ID.

    Args:
            job_id (str): The job ID to fetch the results for.

    Returns:
            Response: The results of the retrosynthesis job.
    """
    if job_id not in results:
        return jsonify({"error": "job not found"}), 404

    if results[job_id]["status"] != "done":
        return jsonify({"error": "job not finished"}), 404

    return jsonify(results.pop(job_id)), 200
