class AiZynthArgs:
	"""Used to create an aizynthfinder instance"""
	def __init__(self, smiles, policy, stock):
		self.smiles = smiles
		# self.config = config_file
		self.policy = []
		self.filter = []
		self.stocks = []
		self.output = None
		self.log_to_file = False
		self.cluster = True
		self.route_distance_model = None
		self.post_processing = True


class RetroRoute:
	"""For a route dictionary, extracts the reaction routes"""
	def __init__(self, route_dict):
		self.reactions = []
		self.depth = 0
		self.previous_depth = 0
		self.route_dict = route_dict
		self.parent = None
		self.parent_smiles = None
		self.node_id_ls = []

	def get_child_smiles(self, node):
		child_smiles_ls = []
		reaction_class_ls = []
		if 'children' in node.keys():
			for child in node['children']:
				if child['type'] == 'mol':
					child_smiles_ls.append(child['smiles'])
				# elif child['type'] == 'reaction':
				else:
					if child['type'] == 'reaction':
						reaction_class = child['metadata']['classification']
						print(reaction_class)
						reaction_class_ls.append(reaction_class)
					if 'children' in child.keys():
						for grandchild in child['children']:
							if grandchild['type'] == 'mol':
								child_smiles_ls.append(grandchild['smiles'])
		return child_smiles_ls, reaction_class_ls

	def find_child_nodes2(self, node):
		if type(node) is dict:
			if node['type'] == 'mol':
				# if dictionary then append to reactions
				child_smiles, reaction_classes = self.get_child_smiles(node)
				reaction_classes = ['' if rxn_class == 'Unassigned' or rxn_class == '0.0 Unrecognized' else rxn_class for rxn_class in reaction_classes]
				self.reactions.append(
					{'node_id': self.node_id(node), 'smiles': node['smiles'], 'node': node, 'depth': self.depth,
					 'parent': self.parent, 'parent_smiles': self.parent_smiles,
					 'child_smiles': child_smiles, 'reaction_class': reaction_classes})
				if self.depth > self.previous_depth:
					self.parent = node
					self.parent_smiles = node['smiles']
			if 'children' in node.keys():
				self.find_child_nodes2(node['children'])
		elif type(node) is list:
			self.previous_depth = self.depth
			self.depth += 1
			for children in node:
				self.find_child_nodes2(children)

	def node_id(self, node):
		"""Creates an id for the node"""
		if self.node_id_ls:
			id_ls = [x.split('node-')[-1] for x in self.node_id_ls]
			# 0-1
			current_depth_id_ls = [int(x.split(f'{self.depth}-')[-1]) for x in id_ls if f'{self.depth}-' in x]
			if current_depth_id_ls:
				next_number = max(current_depth_id_ls) + 1
			else:
				next_number = 0
			# e.g., node-0, node-1-0
			node_id = f'node-{self.depth}-{next_number}'
		else:
			# for first node
			node_id = 'node-0'
		self.node_id_ls.append(node_id)
		return node_id

	def type_check(self, variable):
		if type(variable) is dict:
			return True
		if type(variable) is list:
			return False
