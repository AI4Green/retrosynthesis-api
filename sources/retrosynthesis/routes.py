import os
import json
import time
from flask import request

import sources.retrosynthesis.startup
from sources import app

ACCESS_KEY = os.getenv("KEY", 'retro_key')

@app.route('/', methods=['GET'])
def service_check():
	page_data = {'Message': 'Retrosynthesis service is running', 'Timestamp': time.time()}
	json_dump = json.dumps(page_data)
	return json_dump


@app.route('/retrosynthesis_api/', methods=['GET'])
def retrosynthesis():
	access_key = str(request.args.get('key'))
	if access_key != ACCESS_KEY:
		print("invalid key")
		return json.dumps({'Message': 'invalid key', 'Timestamp': time.time()})
	smiles = str(request.args.get('smiles'))
	solved_route_dict_ls, raw_routes = retrosynthesis_process(smiles)
	page_data = {'Message': solved_route_dict_ls, 'Raw_Routes': raw_routes, 'Timestamp': time.time()}
	json_dump = json.dumps(page_data)
	return json_dump


def retrosynthesis_process(smiles):
	"""
	Takes a smiles string and returns a list of retrosynthetic routes stored as dictionaries
	"""
	# load config containing policy file locations
	print(smiles)
	from rdkit import Chem
	from sources.retrosynthesis.aizynthfinder.interfaces import aizynthcli
	from sources.retrosynthesis.classes import RetroRoute
	# from sources.retrosynthesis.startup import finder
	finder = sources.retrosynthesis.startup.make_config()
	mol = Chem.MolFromSmiles(smiles)
	print(mol)
	aizynthcli._process_single_smiles(smiles, finder, None, False, None, [])
	# Find solved routes and process routes objects into list of dictionaries
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
		route_dic = {'score': route['all_score']['state score'], 'steps': retro_route.reactions,
		             'depth': route['node'].state.max_transforms}
		solved_route_dict.update({f'Route {idx}': route_dic})
	route_dicts = routes.dicts[0:10]
	raw_routes = []
	for idx, route_dict in enumerate(route_dicts, 1):
		raw_routes.append(route_dict)

	return solved_route_dict, raw_routes
