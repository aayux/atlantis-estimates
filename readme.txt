Atlantis: Lattice-based Anonymous Tokens with a Private Metadata Bit

- Security and size estimates derived from [LNP22]
- Base code: https://gitlab.com/ElenaKirshanova/onemoresis_estimates

- Requirements
    * Python 3

- Description of files
Short description of the content:
    * estimator_ntru.py: the NTRU security estimator borrowed from "NTRU Fatigue: How Stretched is Overstretched?", Léo Ducas and Wessel van Woerden.
    * MLWE_security.py, model_BKZ.py, MSIS_security.py, proba.py are MSIS/MLWE security estimator scipts borrowed from "CRYSTALS-Dilithium – Algorithm Specifications and Supporting Documentation".

    The following files have code for estimating signature size and the concrete bit securities
    * cli_que_estim.py: client query size and security estimator
    * iss_tok_estim.py: issuer pretoken size and security estimator 
    * cliq_obt_estim.py: final token size and security estimator

- Experiments
To generate the various estimates run

```
python <experiment-name>_estim.py
```
