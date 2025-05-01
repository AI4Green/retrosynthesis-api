import os
import json
import time
from flask import request

import sources.retrosynthesis.startup
from sources import app

ACCESS_KEY = os.getenv("KEY", 'retro_key')

@app.route('/', methods=['GET'])
def service_check():
	finder = sources.retrosynthesis.startup.make_config()
	# stock1 = finder.config.stock
	page_data = {'Message': 'Retrosynthesis service is running', 'Timestamp': time.time(), 'Version': 'v2.05.25'}
	json_dump = json.dumps(page_data)
	return json_dump


@app.route('/retrosynthesis_api/', methods=['GET'])
def retrosynthesis():
	access_key = str(request.args.get('key'))
	if access_key != ACCESS_KEY:
		print("Invalid key")
		return json.dumps({'Message': 'Invalid key', 'Timestamp': time.time()})

	smiles = str(request.args.get('smiles'))
	enhancement = str(request.args.get('enhancement', 'Default'))
	finder = sources.retrosynthesis.startup.make_config()
	finder.config.search.algorithm_config["enhancement"] = enhancement
	finder.config.search.iteration_limit = int(request.args.get('iterations'))
	print('iterations is:', finder.config.search.iteration_limit)
	finder.config.search.max_transforms = int(request.args.get('max_depth'))
	print('transforms is:', finder.config.search.max_transforms)
	finder.config.search.time_limit = int(request.args.get('time_limit'))
	print('time limit is:', finder.config.search.time_limit)
	solved_route_dict_ls, raw_routes = retrosynthesis_process(smiles, finder)
	print(len(solved_route_dict_ls))
	page_data = {'Message': solved_route_dict_ls, 'Raw_Routes': raw_routes, 'Timestamp': time.time()}

	json_dump = json.dumps(page_data)
	return json_dump


def retrosynthesis_process(smiles, finder):
	"""
    Takes a SMILES string and a pre-configured finder object and returns a list of retrosynthetic routes as dictionaries.
    """
	print(f"Running retrosynthesis for SMILES: {smiles}")

	from rdkit import Chem
	from aizynthfinder.interfaces import aizynthcli
	from sources.retrosynthesis.classes import RetroRoute

	mol = Chem.MolFromSmiles(smiles)
	if not mol:
		raise ValueError("Invalid SMILES string")
	print(f"Molecule generated: {mol}")
	aizynthcli._process_single_smiles(smiles, finder, None, False, None, [], None)
	routes = finder.routes
	solved_routes = []
	for idx, node in enumerate(routes.nodes):
		if node.is_solved is True:
			solved_routes.append(routes[idx])
	solved_routes = solved_routes[0:10]
	solved_route_dict = {}
	for idx, route in enumerate(solved_routes, 1):
		retro_route = RetroRoute(route['dict'])
		retro_route.find_child_nodes2(retro_route.route_dict)
		route_dic = {
			'score': route['all_score']['state score'],
			'steps': retro_route.reactions,
			'depth': route['node'].state.max_transforms,
		}
		solved_route_dict[f"Route {idx}"] = route_dic
	route_dicts = routes.dicts[0:10]
	raw_routes = [route_dict for route_dict in route_dicts]

	return solved_route_dict, raw_routes