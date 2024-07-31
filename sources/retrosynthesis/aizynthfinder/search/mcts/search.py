""" Module containing a class that holds the tree search
"""
from __future__ import annotations

import json
from typing import TYPE_CHECKING
import csv
import networkx as nx
# import matplotlib.pyplot as plt
import random
from ....aizynthfinder.chem import MoleculeDeserializer, MoleculeSerializer
from ....aizynthfinder.search.mcts.node import MctsNode

if TYPE_CHECKING:
    from aizynthfinder.context.config import Configuration
    from aizynthfinder.utils.type_utils import List, Optional

whole_tree = []
network_scores = []
uct_scores = []
class MctsSearchTree:
    """
    Encapsulation of the search tree.

    :ivar root: the root node
    :ivar config: the configuration of the search tree

    :param config: settings of the tree search algorithm
    :param root_smiles: the root will be set to a node representing this molecule, defaults to None
    """

    def __init__(self, config: Configuration, root_smiles: str = None) -> None:

        self.profiling = {
            "expansion_calls": 0,
            "reactants_generations": 0,
            "iterations": 0,
        }

        if root_smiles:
            self.root: Optional[MctsNode] = MctsNode.create_root(
                smiles=root_smiles, tree=self, config=config
            )
        else:
            self.root = None
        self.config = config
        self._graph: Optional[nx.DiGraph] = None

    # @classmethod
    # def from_json(cls, filename: str, config: Configuration) -> "MctsSearchTree":
    #     """
    #     Create a new search tree by deserialization from a JSON file
    #
    #     :param filename: the path to the JSON node
    #     :param config: the configuration of the search
    #     :return: a deserialized tree
    #     """
    #     tree = MctsSearchTree(config)
    #     with open(filename, "r") as fileobj:
    #         dict_ = json.load(fileobj)
    #     mol_deser = MoleculeDeserializer(dict_["molecules"])
    #     tree.root = MctsNode.from_dict(dict_["tree"], tree, config, mol_deser)
    #     return tree

    def backpropagate(self, from_node: MctsNode, value_estimate: float) -> None:
        """
        Backpropagate the value estimate and update all nodes from a
        given node all the way to the root.

        :param from_node: the end node of the route to update
        :param value_estimate: The score value to backpropagate
        """
        current = from_node
        while current is not self.root:
            parent = current.parent
            # For mypy, parent should never by None unless current is the root
            assert parent is not None
            parent.backpropagate(current, value_estimate)
            current = parent

    def graph(self, recreate: bool = False) -> nx.DiGraph:
        """
        Construct a directed graph object with the nodes as
        vertices and the actions as edges attribute "action".

        :param recreate: if True will construct the graph even though it is cached, defaults to False
        :return: the graph object
        :raises ValueError: if the tree is not defined
        """
        if not self.root:
            raise ValueError("Root of search tree is not defined ")

        if not recreate and self._graph:
            return self._graph

        def add_node(node):
            self._graph.add_edge(node.parent, node, action=node.parent[node]["action"])
            for grandchild in node.children:
                add_node(grandchild)

        self._graph = nx.DiGraph()
        # Always add the root
        self._graph.add_node(self.root)
        for child in self.root.children:
            add_node(child)
        return self._graph

    def nodes(self) -> List[MctsNode]:
        """Return all the nodes in the search tree"""
        return list(self.graph())

    def one_iteration(self) -> bool:
        """
        Perform one iteration of
            1. Selection
            2. Expansion
            3. Rollout
            4. Backpropagation

        :return: if a solution was found
        """
        self.profiling["iterations"] += 1
        print(self.profiling['iterations'])
        leaf = self.select_leaf()
        leaf.expand()
        for i in leaf.state.mols:
            print(i.smiles)
        print(leaf.parent)
        print(leaf.children)
        smiles = leaf.state.extract_smiles_from_str_representation()
        current_depth = leaf.state.max_transforms
        solved_mols = []

        counter = 0
        while (not leaf.is_sim_terminal(smiles, self.config.max_transforms, current_depth, self.config)
               and len(smiles) > 0):
            actions, used_mol = leaf.find_actions(smiles, self.config.expansion_policy)
            reactants = None
            print('initial', smiles)
            while actions:
                action = random.choice(actions)
                reactants = leaf.perform_action(action)
                if not reactants:
                    counter +=1
                    actions.remove(action)
            print(actions)

            if used_mol in smiles:
                smiles.remove(used_mol)

            for tree_mol in reactants:
                smiles.append(tree_mol.smiles)


            smiles, solved = leaf.separate_solved_mols(smiles, self.config)
            solved_mols.append(solved)
            current_depth += 1

        if len(smiles) > 0 and len(solved_mols) > 0:
            score = len(solved_mols)/len(smiles)
        if len(solved_mols) == 0:
            score = 0
        if len(smiles) == 0:
            score = 1
        print(leaf.parent.children)
        leaf.backpropagate(leaf, score)
        return leaf.state.is_solved

    def select_leaf(self) -> MctsNode:
        """
        Traverse the tree selecting the most promising child at
        each step until leaf node returned.

        :return: the leaf node
        :raises ValueError: if the tree is not defined
        """
        if not self.root:
            raise ValueError("Root of search tree is not defined ")

        current = self.root
        while current.is_expanded and not current.state.is_solved:
            promising_child = current.promising_child()

            # If promising_child returns None it means that the node
            # is unexpandable, and hence we should break the loop
            if promising_child:
                current = promising_child
        #         current.position = [1, 1, 'select']
        #
        # if current == self.root:
        #     current.position = [0, 0, 'root']
        return current

    # def serialize(self, filename: str) -> None:
    #     """
    #     Serialize the search tree to a JSON file

    #     :param filename: the path to the JSON file
    #     :raises ValueError: if the tree is not defined
    #     """
    #     if not self.root:
    #         raise ValueError("Root of search tree is not defined ")

    #     mol_ser = MoleculeSerializer()
    #     dict_ = {"tree": self.root.serialize(mol_ser), "molecules": mol_ser.store}
    #     with open(filename, "w") as fileobj:
    #         json.dump(dict_, fileobj, indent=2)

    # def create_my_diagram(self, data):
    #     G = nx.DiGraph()

    #     colors = [
    #         '#E6F7FF', '#CCEFFF', '#B3E6FF', '#99DEFF', '#80D4FF',
    #         '#66CBFF', '#4DC1FF', '#33B8FF', '#1AAEFF', '#0094FF',
    #         '#007ACC', '#0061B3', '#004999', '#003280', '#002166',
    #         '#00194D', '#000F33', '#000A1A', '#000000'
    #     ]

    #     for node_data in data:
    #         depth, child_index, operation, parent_pos, smiles, iteration = node_data
    #         parent_pos = tuple(parent_pos) if parent_pos is not None else None

    #         if parent_pos is not None and parent_pos not in G:
    #             G.add_node(parent_pos)
    #         node_id = (depth, child_index)
    #         G.add_node(node_id, operation=operation, smiles=smiles, iteration=iteration)

    #         if parent_pos is not None:
    #             G.add_edge(parent_pos, node_id)

    #     pos = nx.nx_agraph.graphviz_layout(G, prog='dot', args='-Grankdir=TB -Gnodesep=1 -Granksep=2')

    #     max_iteration = max(data, key=lambda x: x[5])[5]
    #     for iteration in range(1, max_iteration + 1):
    #         iter_nodes = [node for node, attrs in G.nodes(data=True) if attrs.get('iteration') == iteration]
    #         color_index = min(iteration - 1, len(colors) - 1)
    #         nx.draw_networkx_nodes(G, pos, nodelist=iter_nodes, node_color=colors[color_index])

    #     nx.draw_networkx_edges(G, pos)
    #     nx.draw_networkx_labels(G, pos, font_size=8, font_color='black')

    #     name = 'iteration' + str(self.profiling["iterations"])
    #     plt.title(name)
    #     plt.savefig(name)
    #     plt.show()

########################################################################################
#                           Storage of outdated functions                              #
########################################################################################

    def one_iteration_od(self) -> bool:

        """
        Perform one iteration of
            1. Selection
            2. Expansion
            3. Rollout
            4. Backpropagation

        :return: if a solution was found
        """
        self.profiling["iterations"] += 1
        print(self.profiling["iterations"])
        leaf = self.select_leaf()
        #leaf.find_position('select', self.profiling["iterations"])
        leaf.expand()
        #whole_tree.append(leaf.position)


        # for value in leaf._children_values:
        #     network_scores.append(value)
        # scores = leaf._children_q() + leaf._children_u()
        # for score in scores:
        #     uct_scores.append(score)

        rollout_child = None
        while not leaf.is_terminal():
            child = leaf.promising_child()
            if not rollout_child:
                rollout_child = child
            if child:
                child.expand()
                #child.find_position('simulate', self.profiling["iterations"])

                #whole_tree.append(child.position)
                #print(whole_tree)

                # for value in leaf._children_values:
                #     network_scores.append(value)
                # scores = leaf._children_q() + leaf._children_u()
                # for score in scores:
                #     uct_scores.append(score)

                leaf = child
        self.backpropagate(leaf, leaf.state.score)

        # it_array = [1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100]
        #
        # if self.profiling["iterations"] in it_array:
        #
        #     self.create_my_diagram(whole_tree)
            # counter_uct = 0
            # total_uct = 0
            # counter_net = 0
            # total_net = 0
            # for score in uct_scores:
            #     if score > 0:
            #         counter_uct +=1
            #         total_uct += score
            # for score in network_scores:
            #     if score > 0:
            #         counter_net +=1
            #         total_net += score

            # print(whole_tree)
            # print(total_uct/counter_uct)
            # print(total_net/counter_net)
            # print(sum(uct_scores)/len(uct_scores))
            # print(sum(network_scores)/len(network_scores))
        return leaf.state.is_solved

# if counter < 1:
#                     with open('policy_testing.csv', 'a') as t:
#                         w = csv.writer(t, lineterminator='\n')                            First iteration policy view
#                         z = child.children_view()
#                         w.writerow(z["actions"])
#                         w.writerow(z["values"])
#                     counter +=1



    # def plot_tree_OD(self, position_list, save_path=None):
    #     G = nx.Graph()

    #     # Adjust the scale factor to control the overall size of the graph
    #     scale_factor = 200

    #     for i, node in enumerate(position_list):
    #         if node:
    #             depth, child_index, _ = node
    #             G.add_node((depth, child_index), pos=(child_index * scale_factor, -depth * scale_factor))

    #     for i, node in enumerate(position_list):
    #         if node:
    #             depth, child_index, _ = node
    #             parent_index = (depth - 1, i - depth + 1)

    #             # Ensure that the parent node exists and has a position assigned
    #             if parent_index in G.nodes and G.nodes[parent_index]:
    #                 G.add_edge(parent_index, (depth, child_index))

    #     pos = nx.get_node_attributes(G, 'pos')

    #     # Ensure that all nodes in the graph have positions assigned
    #     for node in G.nodes:
    #         if node not in pos:
    #             pos[node] = (0, 0)  # Assign a default position for nodes without positions

    #     nx.draw(G, pos, with_labels=True, node_size=700, node_color='skyblue', font_size=8, font_color='black')

    #     if save_path:
    #         plt.savefig(save_path, format='png')

    #     plt.show()


    # def create_tree_diagram_OD(self, data):
    #         # Create an empty graph
    #         G = nx.DiGraph()
    #
    #         # Create a dictionary to keep track of child nodes for each parent position
    #         child_nodes = {}
    #
    #         # Iterate over the nodes in the data
    #         for node_data in data:
    #             depth, child_index, operation, parent_pos, smiles = node_data
    #
    #             # Convert parent position to tuple if not None
    #             parent_pos = tuple(parent_pos) if parent_pos is not None else None
    #
    #         # Check if the parent has already been added to the graph
    #         if parent_pos is not None and parent_pos not in G:
    #             # If not, add the parent node to the graph
    #             G.add_node(parent_pos)
    #
    #         # Add the current node to the graph
    #         node_id = (depth, child_index)
    #         G.add_node(node_id, operation=operation)
    #
    #         # Connect the current node to its parent, if it exists
    #         if parent_pos is not None:
    #             G.add_edge(parent_pos, node_id)
    #
    #         # Check if the child position is already connected to another parent
    #         if parent_pos is not None and node_id in G:
    #             # If the smile strings are different, treat them as separate nodes
    #             if smiles != G.nodes[node_id].get('smiles', ''):
    #                 new_node_id = (depth, child_index, len(G) + 1)  # Use a unique node ID
    #                 G.add_node(new_node_id, operation=operation, smiles=smiles)
    #                 G.add_edge(parent_pos, new_node_id)
    #
    #     # Position the nodes
    #     pos = nx.spring_layout(G, seed=42, k=0.5)
    #
    #     # Draw the graph
    #     nx.draw(G, pos, with_labels=True, node_size=700, node_color='skyblue', font_size=8, font_color='black',
    #             arrows=True)
    #
    #     plt.show()

