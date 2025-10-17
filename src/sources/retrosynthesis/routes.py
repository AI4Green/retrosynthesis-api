from multiprocessing import Queue
import os
import json
import time
from typing import Optional
import uuid
from flask import request, jsonify

import sources.retrosynthesis.startup
from sources import app

# ACCESS_KEYS can be a comma-separated string of keys per client
# e.g. client1-key,client2-key
ACCESS_KEYS = os.getenv("KEYS", "retro_key").split(",")

# These will be reset in app.py to what they need to be
queue: Optional[Queue] = None
results: Optional[dict] = None


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
        return jsonify({"status": "error", "error": "invalid key"}), 403

    # Get job parameters
    smiles = str(request.args.get("smiles"))
    enhancement = str(request.args.get("enhancement", "Default"))
    iteration_limit = int(request.args.get("iterations"))
    max_transforms = int(request.args.get("max_depth"))
    time_limit = int(request.args.get("time_limit"))

    # Create unique job ID
    job_id = str(uuid.uuid4())

    # Queue the analysis
    queue.put(
        (job_id, smiles, enhancement, iteration_limit, max_transforms, time_limit)
    )

    # It takes a long time for the job to appear in results, making it look
    # like an error in the UI. So add the entry here to get ahead
    results[job_id] = {"status": "running", "results": {}}

    return jsonify({"job_id": job_id}), 200


@app.route("/results/<job_id>", methods=["GET"])
def get_results(job_id: str):
    """Get the results for a retrosynthesis job with the given ID.

    Args:
        job_id (str): The job ID to fetch the results for.

    Returns:
        Response: The results of the retrosynthesis job.
    """
    access_key = str(request.args.get("key"))
    if access_key not in ACCESS_KEYS:
        return jsonify({"status": "error", "error": "invalid key"}), 403

    if job_id not in results:
        return jsonify({"status": "error", "error": "job not found"}), 404

    if results[job_id]["status"] == "error":
        # return job information by removing from results
        return jsonify(results.pop(job_id)), 500

    if results[job_id]["status"] == "running":
        # return job information without removing
        return jsonify(results[job_id]), 200

    # return job information by removing from results
    return jsonify(results.pop(job_id)), 200
