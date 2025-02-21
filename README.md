# Retrosynthesis API
A retrosynthesis flask app which acts as an API for AI4Green.

The API uses code and models extracted from the AiZynthFinder package to perform retrosynthesis on a target molecule to predict retrosynthetic routes to the target molecule. - https://github.com/MolecularAI/aizynthfinder

This repository includes a template file and model file. These can be swapped out for alternatives provided that the paths in startup.py are edited to reflect any changes.

**A stock file must be added to sources/retrosynthesis/config_files**.

The API accepts a SMILES string as the target molecule.

The API returns a response with retrosynthetic routes processed for use in AI4Green, unprocessed routes, and an object which can be used to visualise the Monte Carlo Tree Search.

Includes a docker file.

The docker image can be found at: https://hub.docker.com/repository/docker/ai4greeneln/retrosynthesis/general

Supplementary information from our retrosynthesis experiments can be found here: https://doi.org/10.6084/m9.figshare.28450787.v1
