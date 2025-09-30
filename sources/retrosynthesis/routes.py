import os
import json
import time
from flask import request, jsonify
from pathlib import Path

import sources.retrosynthesis.startup
from sources import app

ACCESS_KEY = os.getenv("KEY", "retro_key")
BATCH_OUT_DIR = Path(os.getenv("results", "."))
BATCH_OUT_FILE = BATCH_OUT_DIR / "trees_batch.json"
BATCH_PROGRESS = {}

@app.route("/", methods=["GET"])
def service_check():
    _ = sources.retrosynthesis.startup.make_config()
    page_data = {
        "Message": "Retrosynthesis service is running",
        "Timestamp": time.time(),
        "Version": "v2.05.25",
    }
    return jsonify(page_data)


@app.route("/retrosynthesis_api/", methods=["GET"])
def retrosynthesis():
    access_key = str(request.args.get("key", ""))
    if access_key != ACCESS_KEY:
        return jsonify({"Message": "Invalid key", "Timestamp": time.time()}), 401

    smiles = str(request.args.get("smiles", "")).strip()
    if not smiles:
        return jsonify({"Message": "Missing 'smiles'", "Timestamp": time.time()}), 400

    enhancement = str(request.args.get("enhancement", "Default"))
    iterations = int(request.args.get("iterations", 100))
    max_depth  = int(request.args.get("max_depth", 7))
    time_limit = int(request.args.get("time_limit", 60))

    finder = sources.retrosynthesis.startup.make_config()
    finder.config.search.algorithm_config["enhancement"] = enhancement
    finder.config.search.iteration_limit = iterations
    finder.config.search.max_transforms  = max_depth
    finder.config.search.time_limit      = time_limit

    solved_route_dict_ls, raw_routes = retrosynthesis_process(smiles, finder)
    page_data = {
        "Message": solved_route_dict_ls,
        "Raw_Routes": raw_routes,
        "Timestamp": time.time(),
    }
    return jsonify(page_data)


def retrosynthesis_process(smiles, finder):
    """
    Takes a SMILES string and a pre-configured finder object and returns
    a list of retrosynthetic routes as dictionaries.
    """
    print(f"Running retrosynthesis for SMILES: {smiles}")

    from rdkit import Chem
    from aizynthfinder.interfaces import aizynthcli
    from sources.retrosynthesis.classes import RetroRoute

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError("Invalid SMILES string")

    aizynthcli._process_single_smiles(smiles, finder, None, False, None, [], None)

    routes = finder.routes
    solved_routes = []
    for idx, node in enumerate(routes.nodes):
        if node.is_solved is True:
            solved_routes.append(routes[idx])

    solved_routes = solved_routes[0:10]
    solved_route_dict = {}
    for idx, route in enumerate(solved_routes, 1):
        retro_route = RetroRoute(route["dict"])
        retro_route.find_child_nodes2(retro_route.route_dict)
        route_dic = {
            "score": route["all_score"]["state score"],
            "steps": retro_route.reactions,
            "depth": route["node"].state.max_transforms,
        }
        solved_route_dict[f"Route {idx}"] = route_dic

    route_dicts = routes.dicts[0:10]
    raw_routes = [route_dict for route_dict in route_dicts]

    return solved_route_dict, raw_routes

@app.route("/retrosynthesis_batch_progress/<batch_id>", methods=["GET"])
def retrosynthesis_batch_progress(batch_id):
    info = BATCH_PROGRESS.get(batch_id)
    if not info:
        return jsonify({"status": "unknown"}), 404
    return jsonify(info)

@app.route("/retrosynthesis_batch/", methods=["POST", "OPTIONS"])
def retrosynthesis_batch():
    if request.method == "OPTIONS":
        return ("", 204)

    payload = request.get_json(silent=True) or {}
    access_key = str(payload.get("key", ""))
    if access_key != ACCESS_KEY:
        return jsonify({"Message": "Invalid key", "Timestamp": time.time()}), 401

    smiles_list = payload.get("smiles_list") or []
    if not isinstance(smiles_list, list) or not smiles_list:
        return jsonify({"Message": "smiles_list must be a non-empty list", "Timestamp": time.time()}), 400

    enhancement = str(payload.get("enhancement", "Default"))
    iterations = int(payload.get("iterations", 100))
    max_depth = int(payload.get("max_depth", 7))
    time_limit = int(payload.get("time_limit", 60))

    batch_id = str(payload.get("batch_id", "")) or str(int(time.time()))
    BATCH_PROGRESS[batch_id] = {"total": len(smiles_list), "current": 0, "smiles": "", "status": "running"}

    results = []
    try:
        for idx, smi in enumerate(smiles_list, start=1):
            BATCH_PROGRESS[batch_id].update({"current": idx, "smiles": smi, "status": "running"})

            try:
                finder = sources.retrosynthesis.startup.make_config()
                finder.config.search.algorithm_config["enhancement"] = enhancement
                finder.config.search.iteration_limit = iterations
                finder.config.search.max_transforms = max_depth
                finder.config.search.time_limit = time_limit

                solved_route_dict_ls, raw_routes = retrosynthesis_process(smi, finder)
                results.append(
                    {"smiles": smi, "status": "ok", "routes": solved_route_dict_ls, "raw_routes": raw_routes})
            except Exception as e:
                results.append({"smiles": smi, "status": "failed", "error": str(e), "routes": {}, "raw_routes": []})

        BATCH_PROGRESS[batch_id].update({"status": "done"})
    except Exception:
        BATCH_PROGRESS[batch_id].update({"status": "failed"})
        raise

    batch_entries = []
    for item in results:
        batch_entries.append({
            "target": item["smiles"],
            "status": item.get("status", "ok"),
            "routes": item.get("routes", {}),
            "raw_routes": item.get("raw_routes", []),
            "timestamp": time.time(),
        })
    BATCH_OUT_DIR.mkdir(parents=True, exist_ok=True)
    existing = []
    if BATCH_OUT_FILE.exists():
        try:
            with BATCH_OUT_FILE.open("r", encoding="utf-8") as f:
                data = json.load(f)
                if isinstance(data, dict) and "entries" in data and isinstance(data["entries"], list):
                    existing = data["entries"]
                elif isinstance(data, list):
                    existing = data
        except Exception:
            existing = []
    combined = existing + batch_entries
    tmp_path = BATCH_OUT_FILE.with_suffix(".json.tmp")
    with tmp_path.open("w", encoding="utf-8") as f:
        json.dump(combined, f, ensure_ascii=False, indent=2)
    tmp_path.replace(BATCH_OUT_FILE)

    return jsonify(
        {"Message": "Batch completed", "results": results, "Timestamp": time.time()}
    )
