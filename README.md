### Examples of SMEFTatNLO computations

Standalone bash scripts are the following:

- `gen_single_top_4f.sh`: Compute single top production at NLO in QCD with `cQq13` or `cQq83` four-fermion operator (at fixed order, in the five flavour scheme where t and s channels are not separated). NLO overlap with tW associated production is removed with the `$$` syntax. The W width has to be set to zero to preserve gauge invariance (and ensure pole cancellation). In the `cQq83` case, a user loop filter is needed to include loops involving both a gluon and a W that are otherwise not considered as pure QCD ones.

- `test-topfcnc-pp2th.sh`: Compute FCNC production of a top and a Higgs (`pp>th`) at NLO in QCD with `RCtphi` operator coefficient from the [TopFCNC](https://feynrules.irmp.ucl.ac.be/wiki/TopFCNC) model. The `$$ t~` syntax is used to remove the `pp>tt~` contribution with `t~>hj` which arises at NLO.

- `gen-gen-pp2tta-pta.sh`: Compute top-pair production in association with a photon (`pp>tta`) with `ctZ` dipole dependence. Produce `pT(γ)` distributions for the SM, linear, and quadratic terms separately as well as the associated differential k-factor plots, using fixed-order mode.
