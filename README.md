# Retrosynthesis API
A retrosynthesis flask app which acts as an API for AI4Green.

The API uses AiZynthFinder to perform retrosynthesis on a target molecule to predict retrosynthetic routes to the target molecule. - https://github.com/MolecularAI/aizynthfinder

The API accepts a SMILES string as the target molecule.

The API returns a response with retrosynthetic routes processed for use in AI4Green, unprocessed routes, and an object which can be used to visualise the Monte Carlo Tree Search.

Includes a docker file.

The docker image can be found at: https://hub.docker.com/repository/docker/ai4greeneln/retrosynthesis/general
