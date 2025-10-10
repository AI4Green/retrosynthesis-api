from multiprocessing import Process, Queue, Manager
import time, uuid, os

_manager = Manager()
results = _manager.dict()
_queue = Queue()


def retrosynthesis_process(smiles, finder):
	"""
    Takes a SMILES string and a pre-configured finder object and returns a list of retrosynthetic routes as dictionaries.
    """

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