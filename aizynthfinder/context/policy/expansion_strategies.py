""" Module containing classes that implements different expansion policy strategies
"""

from __future__ import annotations

import abc
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from aizynthfinder.chem import SmilesBasedRetroReaction, TemplatedRetroReaction
from aizynthfinder.context.policy.utils import _make_fingerprint
from aizynthfinder.utils.exceptions import PolicyException
from aizynthfinder.utils.logging import logger
from aizynthfinder.utils.models import load_model
from aizynthfinder.context.policy.post_sm_grouping import compute_center_group_key


if TYPE_CHECKING:
    from aizynthfinder.chem import TreeMolecule
    from aizynthfinder.chem.reaction import RetroReaction
    from aizynthfinder.context.config import Configuration
    from aizynthfinder.utils.type_utils import (
        Any,
        Dict,
        List,
        Optional,
        Sequence,
        StrDict,
        Tuple,
    )


class ExpansionStrategy(abc.ABC):
    """
    A base class for all expansion strategies.

    The strategy can be used by either calling the `get_actions` method
    of by calling the instantiated class with a list of molecule.

    .. code-block::

        expander = MyExpansionStrategy("dummy", config)
        actions, priors = expander.get_actions(molecules)
        actions, priors = expander(molecules)

    :param key: the key or label
    :param config: the configuration of the tree search
    """

    _required_kwargs: List[str] = []

    def __init__(self, key: str, config: Configuration, **kwargs: str) -> None:
        if any(name not in kwargs for name in self._required_kwargs):
            raise PolicyException(
                f"A {self.__class__.__name__} class needs to be initiated "
                f"with keyword arguments: {', '.join(self._required_kwargs)}"
            )
        self._config = config
        self._logger = logger()
        self.key = key


    def __call__(
        self,
        molecules: Sequence[TreeMolecule],
        cache_molecules: Optional[Sequence[TreeMolecule]] = None,
    ) -> Tuple[List[RetroReaction], List[float]]:
        return self.get_actions(molecules, cache_molecules)

    @abc.abstractmethod
    def get_actions(
        self,
        molecules: Sequence[TreeMolecule],
        cache_molecules: Optional[Sequence[TreeMolecule]] = None,
    ) -> Tuple[List[RetroReaction], List[float]]:
        """
        Get all the probable actions of a set of molecules

        :param molecules: the molecules to consider
        :param cache_molecules: additional molecules to submit to the expansion
                                  policy but that only will be cached for later use
        :return: the actions and the priors of those actions
        """

    def reset_cache(self) -> None:
        """Reset the prediction cache"""


class MultiExpansionStrategy(ExpansionStrategy):
    """
    A base class for combining multiple expansion strategies.

    The strategy can be used by either calling the `get_actions` method
    or by calling the instantiated class with a list of molecules.

    :ivar expansion_strategy_keys: the keys of the selected expansion strategies
    :ivar additive_expansion: a conditional setting to specify whether all the actions
        and priors of the selected expansion strategies should be combined or not.
        Defaults to False.
    :ivar expansion_strategy_weights: a list of weights for each expansion strategy.
        The weights should sum to one. Exception is the default, where unity weight
        is associated to each strategy.

    :param key: the key or label
    :param config: the configuration of the tree search
    :param expansion_strategies: the keys of the selected expansion strategies. All keys
        of the selected expansion strategies must exist in the expansion policies listed
        in config
    """

    _required_kwargs = ["expansion_strategies"]

    def __init__(
        self,
        key: str,
        config: Configuration,
        **kwargs: Any,
    ) -> None:
        super().__init__(key, config, **kwargs)
        self._config = config
        self._expansion_strategies: List[ExpansionStrategy] = []
        self.expansion_strategy_keys = kwargs["expansion_strategies"]

        self.cutoff_number = kwargs.get("cutoff_number")
        if self.cutoff_number:
            print(f"Setting multi-expansion cutoff_number: {self.cutoff_number}")

        self.expansion_strategy_weights = self._set_expansion_strategy_weights(kwargs)
        self.additive_expansion: bool = bool(kwargs.get("additive_expansion", False))
        self._logger.info(
            f"Multi-expansion strategy with policies: {self.expansion_strategy_keys}"
            f", and corresponding weights: {self.expansion_strategy_weights}"
        )

    def get_actions(
        self,
        molecules: Sequence[TreeMolecule],
        cache_molecules: Optional[Sequence[TreeMolecule]] = None,
    ) -> Tuple[List[RetroReaction], List[float]]:
        """
        Get all the probable actions of a set of molecules, using the selected policies.

        The default implementation combines all the actions and priors of the
        selected expansion strategies into two lists respectively if the
        'additive_expansion' setting is set to True. This function can be overridden by
        a sub class to combine different expansion strategies in different ways.

        :param molecules: the molecules to consider
        :param cache_molecules: additional molecules to submit to the expansion
            policy but that only will be cached for later use
        :return: the actions and the priors of those actions
        :raises: PolicyException: if the policy isn't selected
        """
        expansion_strategies = self._get_expansion_strategies_from_config()

        all_possible_actions = []
        all_priors = []
        for expansion_strategy, expansion_strategy_weight in zip(
            expansion_strategies, self.expansion_strategy_weights
        ):
            possible_actions, priors = expansion_strategy.get_actions(
                molecules, cache_molecules
            )

            all_possible_actions.extend(possible_actions)
            if not self.additive_expansion and all_possible_actions:
                all_priors.extend(priors)
                break

            weighted_prior = [expansion_strategy_weight * p for p in priors]

            all_priors.extend(weighted_prior)

        all_possible_actions, all_priors = self._prune_actions(
            all_possible_actions, all_priors
        )
        return all_possible_actions, all_priors

    def _get_expansion_strategies_from_config(self) -> List[ExpansionStrategy]:
        if self._expansion_strategies:
            return self._expansion_strategies

        if not all(
            key in self._config.expansion_policy.items
            for key in self.expansion_strategy_keys
        ):
            raise ValueError(
                "The input expansion strategy keys must exist in the "
                "expansion policies listed in config"
            )
        self._expansion_strategies = [
            self._config.expansion_policy[key] for key in self.expansion_strategy_keys
        ]

        for expansion_strategy, weight in zip(
            self._expansion_strategies, self.expansion_strategy_weights
        ):
            if not getattr(expansion_strategy, "rescale_prior", True) and weight < 1:
                setattr(expansion_strategy, "rescale_prior", True)
                self._logger.info(
                    f"Enforcing {expansion_strategy.key}.rescale_prior=True"
                )
        return self._expansion_strategies

    def _prune_actions(
        self, actions: List[RetroReaction], priors: List[float]
    ) -> Tuple[List[RetroReaction], List[float]]:
        """
        Prune the actions if a maximum number of actions is specified.

        :param actions: list of predicted actions
        :param priors: list of prediction probabilities
        :return: the top 'self.cutoff_number' actions and corresponding priors.
        """
        if not self.cutoff_number:
            return actions, priors

        sortidx = np.argsort(np.array(priors))[::-1].astype(int)
        priors = [priors[idx] for idx in sortidx[0 : self.cutoff_number]]
        actions = [actions[idx] for idx in sortidx[0 : self.cutoff_number]]
        return actions, priors

    def _set_expansion_strategy_weights(self, kwargs: StrDict) -> List[float]:
        """
        Set the weights of each expansion strategy using the input kwargs from config.
        The weights in the config should sum to one.
        If not set in the config file, the weights default to one for each strategy
        (for backwards compatibility).

        :param kwargs: input arguments to the MultiExpansionStrategy
        :raises: ValueError if weights from the config file do not sum to one.
        :return: a list of expansion strategy weights
        """
        if not "expansion_strategy_weights" in kwargs:
            return [1.0 for _ in self.expansion_strategy_keys]

        expansion_strategy_weights = kwargs["expansion_strategy_weights"]
        sum_weights = sum(expansion_strategy_weights)

        if sum_weights != 1:
            raise ValueError(
                "The expansion strategy weights in MultiExpansion should "
                "sum to one. -> "
                f"sum({expansion_strategy_weights})={sum_weights}."
            )

        return expansion_strategy_weights


class TemplateBasedExpansionStrategy(ExpansionStrategy):
    """
    A template-based expansion strategy that will return `TemplatedRetroReaction` objects upon expansion.

    :ivar template_column: the column in the template file that contains the templates
    :ivar cutoff_cumulative: the accumulative probability of the suggested templates
    :ivar cutoff_number: the maximum number of templates to returned
    :ivar use_rdchiral: a boolean to apply templates with RDChiral
    :ivar use_remote_models: a boolean to connect to remote TensorFlow servers
    :ivar rescale_prior: a boolean to apply softmax to the priors
    :ivar chiral_fingerprints: if True will base expansion on chiral fingerprint
    :ivar mask: a boolean vector of masks for the reaction templates. The length of the vector should be equal to the
        number of templates. It is set to None if no mask file is provided as input.

    :param key: the key or label
    :param config: the configuration of the tree search
    :param model: the source of the policy model
    :param template: the path to a HDF5 file with the templates
    :raises PolicyException: if the length of the model output vector is not same as the
        number of templates
    """

    _required_kwargs = [
        "model",
        "template",
    ]

    def __init__(self, key: str, config: Configuration, **kwargs: str) -> None:
        super().__init__(key, config, **kwargs)

        source = kwargs["model"]
        templatefile = kwargs["template"]
        maskfile: str = kwargs.get("mask", "")
        self.template_column: str = kwargs.get("template_column", "retro_template")
        self.cutoff_cumulative: float = float(kwargs.get("cutoff_cumulative", 0.995))
        self.cutoff_number: int = int(kwargs.get("cutoff_number", 50))
        self.use_rdchiral: bool = bool(kwargs.get("use_rdchiral", True))
        self.use_remote_models: bool = bool(kwargs.get("use_remote_models", False))
        self.rescale_prior: bool = bool(kwargs.get("rescale_prior", False))
        self.chiral_fingerprints = bool(kwargs.get("chiral_fingerprints", False))

        self._logger.info(
            f"Loading template-based expansion policy model from {source} to {self.key}"
        )
        self.model = load_model(source, self.key, self.use_remote_models)

        self._logger.info(f"Loading templates from {templatefile} to {self.key}")
        if templatefile.endswith(".csv.gz") or templatefile.endswith(".csv"):
            self.templates: pd.DataFrame = pd.read_csv(
                templatefile, index_col=0, sep="\t"
            )
        else:
            self.templates = pd.read_hdf(templatefile, "table")

        self.mask: Optional[np.ndarray] = (
            self._load_mask_file(maskfile) if maskfile else None
        )

        if hasattr(self.model, "output_size") and len(self.templates) != self.model.output_size:
            raise PolicyException(
                f"The number of templates ({len(self.templates)}) does not agree with the "  
                f"output dimensions of the model ({self.model.output_size})"
            )
        self._cache: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
        self._cutoff_impl_name: str = str(kwargs.get("cutoff_impl", "templates")).lower()
        self._cutoff_impl = getattr(self, f"_cutoff_predictions_{self._cutoff_impl_name}", None)
        if self._cutoff_impl is None:
            raise PolicyException(f"Unknown cutoff_impl='{self._cutoff_impl_name}'. Use 'templates' or 'groups'.")

        self._group_index = None
        if self._cutoff_impl_name == "groups":
            self._build_group_index()

    def get_actions(
        self,
        molecules: Sequence[TreeMolecule],
        cache_molecules: Optional[Sequence[TreeMolecule]] = None,
    ) -> Tuple[List[RetroReaction], List[float]]:
        """
        Get all the probable actions of a set of molecules, using the selected policies and given cutoffs

        :param molecules: the molecules to consider
        :param cache_molecules: additional molecules to submit to the expansion
                                  policy but that only will be cached for later use
        :return: the actions and the priors of those actions
        """

        possible_actions = []
        priors: List[float] = []
        cache_molecules = cache_molecules or []
        self._update_cache(list(molecules) + list(cache_molecules))

        for mol in molecules:
            probable_transforms_idx, probs = self._cache[mol.inchi_key]
            possible_moves = self.templates.iloc[probable_transforms_idx]
            if self.rescale_prior:
                probs /= probs.sum()
            priors.extend(probs)
            for idx, (move_index, move) in enumerate(possible_moves.iterrows()):
                metadata = dict(move)
                del metadata[self.template_column]
                metadata["policy_probability"] = float(probs[idx].round(4))
                metadata["policy_probability_rank"] = idx
                metadata["policy_name"] = self.key
                metadata["template_code"] = move_index
                metadata["template"] = move[self.template_column]
                possible_actions.append(
                    TemplatedRetroReaction(
                        mol,
                        smarts=move[self.template_column],
                        metadata=metadata,
                        use_rdchiral=self.use_rdchiral,
                    )
                )
        return possible_actions, priors  # type: ignore

    def reset_cache(self) -> None:
        """Reset the prediction cache"""
        self._cache = {}

    def _build_group_index(self) -> None:
        """
        Build a mapping: group_key -> list of row positions of N templates.
        Uses the template SMARTS to compute a product-side reacting-center + first shell signature.
        """
        df = self.templates

        if self.template_column in df.columns:
            smarts_col = self.template_column
        else:
            smarts_col = None
            for c in ("retro_template", "template", "smarts", "reaction_smarts"):
                if c in df.columns:
                    smarts_col = c
                    break
            if smarts_col is None:
                self._logger.info("Grouping: no SMARTS column found; skipping")
                self._group_index = None
                return

        print(f"[DEBUG groups] computing group keys for {len(df)} templates "
              f"using column '{smarts_col}'")
        if "group_key" not in df.columns:
            print("[DEBUG groups] no existing group_key column; computing…")
            df["group_key"] = df[smarts_col].fillna("").map(compute_center_group_key)

        from collections import defaultdict
        gi = defaultdict(list)
        for pos, g in enumerate(df["group_key"].values):
            gi[g].append(pos)
        self._group_index = gi

        ng = len(gi)
        avg = len(df) / max(ng, 1)
        print(f"[DEBUG groups] built {ng} groups over {len(df)} templates "
              f"(avg {avg:.2f}/group)")

        sizes = sorted(((k, len(v)) for k, v in gi.items()), key=lambda x: x[1], reverse=True)[:10]
        if sizes:
            show = ", ".join([f"{k[:10]}…:{n}" if isinstance(k, str) else f"{str(k)[:10]}…:{n}" for k, n in sizes])
            print(f"[DEBUG groups] largest groups (key:count): {show}")


    def _dbg_row_info(self, pos: int) -> str:
        """Build a short string with template identity for position `pos`."""
        df = self.templates
        try:
            label = df.index.values[pos]
        except Exception:
            label = pos
        get = lambda col: (str(df.iloc[pos][col]) if col in df.columns else "")
        thash = get("template_hash")
        gkey  = get("group_key")
        # abbreviate long hashes
        if isinstance(thash, str) and len(thash) > 10:
            thash = thash[:10] + "…"
        if isinstance(gkey, str) and len(gkey) > 10:
            gkey = gkey[:10] + "…"
        return f"pos={pos} label={label} hash={thash} group={gkey}"

    def _dbg_dump_selection(self, preds: np.ndarray, selected: np.ndarray, tag: str) -> None:
        """Print the final selected templates (bounded)."""
        if selected is None or len(selected) == 0:
            print(f"[DEBUG cutoff:{tag}] no templates selected")
            return
        topn = min(25, len(selected))
        total = float(np.nansum(preds))
        cum = float(np.nansum(preds[selected]))
        print(f"[DEBUG cutoff:{tag}] selected {len(selected)} templates "
              f"(mass {cum:.6g} of {total:.6g}); showing top {topn}")
        # order by prob descending
        order = np.argsort(preds[selected])[::-1]
        for rank, rel in enumerate(order[:topn], start=1):
            pos = int(selected[rel])
            p = float(preds[pos])
            print(f"  #{rank:>2}  p={p:.6g}  {self._dbg_row_info(pos)}")

    def _cutoff_predictions(self, predictions: np.ndarray) -> np.ndarray:
        """
        Select up to TOTAL_TARGET templates:
          If GROUPING:
            1) Rank groups by total probability; keep top NUM_GROUPS groups.
            2) From each kept group, take PER_GROUP templates
            3) Top up to TOTAL_TARGET with global best
          Else:
            Take global top TOTAL_TARGET templates.
        """

        # ----------------------------------
        GROUPING = True
        NUM_GROUPS = 20
        PER_GROUP = 5
        TOTAL_TARGET = 100
        FILL_FROM_GLOBAL = True
        # ----------------------------------

        preds = predictions.copy().astype(float)

        if getattr(self, "mask", None) is not None:
            preds[~self.mask] = 0.0

        all_idx = np.arange(preds.size, dtype=np.int32)
        global_order = all_idx[np.argsort(preds)[::-1]]
        if not GROUPING:
            return global_order[: min(TOTAL_TARGET, global_order.size)]
        if getattr(self, "_group_index", None) is None:
            self._build_group_index()
        gi = getattr(self, "_group_index", None)

        if not gi or NUM_GROUPS <= 0 or PER_GROUP <= 0:
            return global_order[: min(TOTAL_TARGET, global_order.size)]

        group_items = list(gi.items())
        group_sums = np.array(
            [float(preds[np.asarray(members, dtype=np.int32)].sum()) for _, members in group_items],
            dtype=float,
        )

        # Order groups by total mass
        g_order = np.argsort(group_sums)[::-1]

        selected: list[int] = []

        # Take top NUM_GROUPS groups
        for g_pos in g_order[: min(NUM_GROUPS, len(g_order))]:
            _, members = group_items[g_pos]
            members = np.asarray(members, dtype=np.int32)
            if members.size == 0:
                continue

            # Top PER_GROUP within this group
            local_sorted = members[np.argsort(preds[members])[::-1]]
            take_k = int(min(PER_GROUP, local_sorted.size))

            count = 0
            for idx in local_sorted:
                if count >= take_k:
                    break
                if preds[idx] > 0.0:
                    selected.append(int(idx))
                    count += 1

            if len(selected) >= TOTAL_TARGET:
                break

        # Top up from global templates if needed
        if FILL_FROM_GLOBAL and len(selected) < TOTAL_TARGET:
            chosen = set(selected)
            for idx in global_order:
                if len(selected) >= TOTAL_TARGET:
                    break
                if idx not in chosen and preds[idx] > 0.0:
                    selected.append(int(idx))

        if not selected:
            return global_order[: min(TOTAL_TARGET, global_order.size)]
        if len(selected) > TOTAL_TARGET:
            selected = selected[:TOTAL_TARGET]

        return np.asarray(selected, dtype=np.int32)

    def _cutoff_predictions_templates(self, predictions: np.ndarray) -> np.ndarray:
        """
        Get the top transformations, by selecting those that have:
            * cumulative probability less than a threshold (cutoff_cumulative)
            * or at most N (cutoff_number)
        """
        # # DEBUG
        # TARGET_HASH = "0c33fbaaea663fe5aa40102c16fffe484773b331031394dbfad9e1eace47a3c2"
        #
        # df = None
        # if hasattr(self, "templates"):
        #     df = self.templates
        # elif hasattr(self, "_template_library") and hasattr(self._template_library, "data"):
        #     df = self._template_library.data
        #
        # if df is not None:
        #     row = df.loc[df["template_hash"] == TARGET_HASH]
        #     if len(row) == 1:
        #         try:
        #             iloc = int(df.index.get_indexer([row.index[0]])[0])
        #         except Exception:
        #             code = int(row["template_code"].iloc[0])
        #             iloc = int(df.index.get_indexer([code])[0])
        #         p = float(predictions[iloc])
        #         rank = int((predictions > p).sum() + 1)
        #         pct = 100.0 * rank / predictions.size
        #         masked = (getattr(self, "mask", None) is not None and not self.mask[iloc])
        #         print(f"[POLICY] prob={p:.6g} rank={rank}/{predictions.size} (top {pct:.4f}%)"
        #               + (" [MASKED]" if masked else ""))
        # # END

        preds = predictions.copy()
        if self.mask is not None:
            preds[~self.mask] = 0

        sortidx = np.argsort(preds)[::-1]
        cumsum: np.ndarray = np.cumsum(preds[sortidx])

        if any(cumsum >= self.cutoff_cumulative):
            maxidx = int(np.argmin(cumsum < self.cutoff_cumulative))
        else:
            maxidx = len(cumsum)

        maxidx = min(maxidx, self.cutoff_number) or 1

        selected = sortidx[:maxidx]
        try:
            self._dbg_dump_selection(preds, selected, tag="templates")
        except Exception as _e:
            print(f"[DEBUG cutoff:templates] printing failed: {_e}")
        return selected

    def _cutoff_predictions_groups(self, predictions: np.ndarray) -> np.ndarray:
        """
        Group-level cutoff:
          1) Sum probabilities per group.
          2) Select top groups by cumulative probability and cutoff_number (applied to groups).
          3) Return template indices from selected groups, sorted by their individual prob,
             then truncate to self.cutoff_number templates.
        """

        if getattr(self, "_group_index", None) is None:
            self._build_group_index()

        preds = predictions.copy()
        if self.mask is not None:
            preds[~self.mask] = 0

        gi = getattr(self, "_group_index", None)
        if not gi:
            print("[DEBUG cutoff:groups] no group index; falling back to template cutoff")
            return self._cutoff_predictions_templates(preds)

        group_keys = list(gi.keys())
        group_members = [gi[k] for k in group_keys]
        group_sums = np.array([float(preds[m].sum()) for m in group_members], dtype=float)

        g_order = np.argsort(group_sums)[::-1]
        g_cumsum = np.cumsum(group_sums[g_order])

        if any(g_cumsum >= self.cutoff_cumulative):
            gmax = int(np.argmin(g_cumsum < self.cutoff_cumulative))
        else:
            gmax = len(g_order)

        gmax = min(gmax, self.cutoff_number) or 1

        print(f"[DEBUG cutoff:groups] considering {len(group_keys)} groups; "
              f"selecting top {gmax} by group mass (cum ≤ {self.cutoff_cumulative})")
        show_g = min(10, len(g_order))
        for rank, gi_pos in enumerate(g_order[:show_g], start=1):
            key = group_keys[gi_pos]
            members = group_members[gi_pos]
            mass = float(group_sums[gi_pos])
            print(f"  [group #{rank:>2}] mass={mass:.6g}  key={str(key)[:10]}…  |members|={len(members)}")

        kept_indices = []
        for i in g_order[:gmax]:
            kept_indices.extend(group_members[i])

        if not kept_indices:
            print("[DEBUG cutoff:groups] kept_indices empty; falling back to template cutoff")
            return self._cutoff_predictions_templates(preds)

        kept_indices = np.array(kept_indices, dtype=np.int32)
        order = np.argsort(preds[kept_indices])[::-1]
        out = kept_indices[order]

        if len(out) > self.cutoff_number:
            out = out[: self.cutoff_number]

        try:
            self._dbg_dump_selection(preds, out, tag="groups")
        except Exception as _e:
            print(f"[DEBUG cutoff:groups] printing failed: {_e}")

        return out

    def _update_cache(self, molecules: Sequence[TreeMolecule]) -> None:
        pred_inchis = []
        fp_list = []
        for molecule in molecules:
            if molecule.inchi_key in self._cache or molecule.inchi_key in pred_inchis:
                continue
            fp_list.append(
                _make_fingerprint(molecule, self.model, self.chiral_fingerprints)
            )
            pred_inchis.append(molecule.inchi_key)

        if not pred_inchis:
            return

        pred_list = np.asarray(self.model.predict(np.vstack(fp_list)))
        for pred, inchi in zip(pred_list, pred_inchis):
            probable_transforms_idx = self._cutoff_predictions(pred)
            self._cache[inchi] = (
                probable_transforms_idx,
                pred[probable_transforms_idx],
            )


class TemplateBasedDirectExpansionStrategy(TemplateBasedExpansionStrategy):
    """
    A template-based expansion strategy that will return `SmilesBasedRetroReaction` objects upon expansion
    by directly applying the template

    :param key: the key or label
    :param config: the configuration of the tree search
    :param source: the source of the policy model
    :param templatefile: the path to a HDF5 file with the templates
    :raises PolicyException: if the length of the model output vector is not same as the number of templates
    """

    def get_actions(
        self,
        molecules: Sequence[TreeMolecule],
        cache_molecules: Optional[Sequence[TreeMolecule]] = None,
    ) -> Tuple[List[RetroReaction], List[float]]:
        """
        Get all the probable actions of a set of molecules, using the selected policies and given cutoffs

        :param molecules: the molecules to consider
        :param cache_molecules: additional molecules to submit to the expansion
            policy but that only will be cached for later use
        :return: the actions and the priors of those actions
        """
        possible_actions = []
        priors = []

        super_actions, super_priors = super().get_actions(molecules, cache_molecules)
        for templated_action, prior in zip(super_actions, super_priors):
            for reactants in templated_action.reactants:
                reactants_str = ".".join(mol.smiles for mol in reactants)
                new_action = SmilesBasedRetroReaction(
                    templated_action.mol,
                    metadata=templated_action.metadata,
                    reactants_str=reactants_str,
                )
                possible_actions.append(new_action)
                priors.append(prior)

        return possible_actions, priors  # type: ignore
