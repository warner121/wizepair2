# Matched Molecular Pairs Analysis (MMPA) based upon Maximum Common Substructure (MCS)
Python (re-)implementation of the algorithm described in: **WizePairZ: A Novel Algorithm to Identify, Encode, and Exploit Matched Molecular Pairs with Unspecified Cores in Medicinal Chemistry** by Daniel J. Warner, Edward J. Griffen, and Stephen A. St-Gallay: *J. Chem. Inf. Model. 2010, 50, 8, 1350–1357, August 6, 2010.*

### Notebooks

The best way to understand how this library can be used is to look at the [examples](https://github.com/warner121/wizepair2/blob/master/notebooks/examples.ipynb).

There are a number of other notebooks illustrating more complex workflows, including a complete reproduction of the [histone deacetylase 1 example](https://github.com/warner121/wizepair2/tree/master/notebooks/hdac1) that was described in the [original publication](https://pubs.acs.org/doi/10.1021/ci100084s).

### Installation

This implementation makes extensive use of the open source cheminformatics toolkit [RDKit](https://rdkit.org/) and the graph toolkit [networkx](https://networkx.org/). Both packages can be installed, along with all other dependencies, by poetry as below.

1. Install poetry (if you don't already have it):
```shell
curl -sSL https://install.python-poetry.org | python3 -
```
2. Clone the wisepair2 repository:  
```shell
git clone https://github.com/warner121/wizepair2.git
cd wizepair2
```
3. Install dependencies and the project in editable mode:
```shell
poetry install
```
4. Run the unit tests.
```shell
poetry run pytest tests/wizepair2/test_mmp.py
```
```text
======================================================== test session starts =========================================================
platform darwin -- Python 3.12.5, pytest-8.3.5, pluggy-1.5.0
rootdir: /Users/daniel/wizepair2
configfile: pyproject.toml
plugins: anyio-4.9.0
collected 14 items                                                                                                                   

tests/wizepair2/test_mmp.py ..............                                                                                     [100%]

======================================================== 14 passed in 37.39s =========================================================
```

### Starting the server

The core functionality is exposed through a fastapi webservice, so that it may be hosted scalable cloud infrastructure. To start the service on your local machine, simply type:

```shell
poetry run python main 
```

Alternatively, you can build the docker image and run the service in a container:
```shell
docker build -t wizepair2 .
docker run -p 8000:8000 wizepair2
```

### Execution

SMILES strings must be provided in the request body as an array of dictionaries, as illustrated below. There is no limit (other than memory and processing time) to the number or length of the strings you pass. In another terminal:

```shell
curl --location --request POST 'http://127.0.0.1:8000/wizepair2/api/v1.0/mmp?strictness=7' \
--header 'Content-Type: application/json' \
--data-raw '[
    {
        "smiles1": "c1ccccc1",
        "smiles2": "c1ccccn1"
    },
    {
        "smiles1": "c1ccccc1",
        "smiles2": "c1cnccn1"
    }
]'
```

The request above should yield the following response. This contains 4 attempted SMIRKS encodings for each pair of molecules, with decreasing atomic environment radii as described in the original paper. The response also contains sanitised fragments describing the transformation, the percent common structure between the two molecules, and a boolean indicating if the SMIRKS was valid. Validity is established when application of the `smirks` to `smiles1` successfully yields `smiles2` (potentially amongst other products).

```json
[
    [
        {
            "fragment1": "[cH2][cH][cH][cH][cH2]",
            "fragment2": "[cH2][cH][n][cH][cH2]",
            "percentmcs": 0.8333333333333334,
            "radius": 4,
            "smiles1": "c1ccccc1",
            "smiles2": "c1ccncc1",
            "smirks": "[#6:5](:[#6:4](:[#6:3](:[#6:2](:[#6:1]-[H])-[H])-[H])-[H])-[H]>>[#6:5](:[#6:4](:[#7:3]:[#6:2](:[#6:1]-[H])-[H])-[H])-[H]",
            "valid": true
        },
        {
            "fragment1": "[cH2][cH][cH][cH][cH2]",
            "fragment2": "[cH2][cH][n][cH][cH2]",
            "percentmcs": 0.8333333333333334,
            "radius": 3,
            "smiles1": "c1ccccc1",
            "smiles2": "c1ccncc1",
            "smirks": "[#6:5](:[#6:4](:[#6:3](:[#6:2](:[#6:1]-[H])-[H])-[H])-[H])-[H]>>[#6:5](:[#6:4](:[#7:3]:[#6:2](:[#6:1]-[H])-[H])-[H])-[H]",
            "valid": true
        },
        {
            "fragment1": "[cH2][cH][cH][cH][cH2]",
            "fragment2": "[cH2][cH][n][cH][cH2]",
            "percentmcs": 0.8333333333333334,
            "radius": 2,
            "smiles1": "c1ccccc1",
            "smiles2": "c1ccncc1",
            "smirks": "[#6:5](:[#6:4](:[#6:3](:[#6:2](:[#6:1]-[H])-[H])-[H])-[H])-[H]>>[#6:5](:[#6:4](:[#7:3]:[#6:2](:[#6:1]-[H])-[H])-[H])-[H]",
            "valid": true
        },
        {
            "fragment1": "[cH2][cH][cH2]",
            "fragment2": "[cH2][n][cH2]",
            "percentmcs": 0.8333333333333334,
            "radius": 1,
            "smiles1": "c1ccccc1",
            "smiles2": "c1ccncc1",
            "smirks": "[#6:4](:[#6:3](:[#6:2]-[H])-[H])-[H]>>[#6:4](:[#7:3]:[#6:2]-[H])-[H]",
            "valid": true
        }
    ],
    [
        {
            "fragment1": "[cH]1[cH][cH][cH][cH][cH]1",
            "fragment2": "[cH]1[cH][n][cH][cH][n]1",
            "percentmcs": 0.6666666666666666,
            "radius": 4,
            "smiles1": "c1ccccc1",
            "smiles2": "c1cnccn1",
            "smirks": "[#6:6]1(:[#6:5](:[#6:4](:[#6:3](:[#6:2](:[#6:1]:1-[H])-[H])-[H])-[H])-[H])-[H]>>[#6:6]1(:[#6:5](:[#7:4]:[#6:3](:[#6:2](:[#7:1]:1)-[H])-[H])-[H])-[H]",
            "valid": true
        },
        {
            "fragment1": "[cH]1[cH][cH][cH][cH][cH]1",
            "fragment2": "[cH]1[cH][n][cH][cH][n]1",
            "percentmcs": 0.6666666666666666,
            "radius": 3,
            "smiles1": "c1ccccc1",
            "smiles2": "c1cnccn1",
            "smirks": "[#6:6]1(:[#6:5](:[#6:4](:[#6:3](:[#6:2](:[#6:1]:1-[H])-[H])-[H])-[H])-[H])-[H]>>[#6:6]1(:[#6:5](:[#7:4]:[#6:3](:[#6:2](:[#7:1]:1)-[H])-[H])-[H])-[H]",
            "valid": true
        },
        {
            "fragment1": "[cH]1[cH][cH][cH][cH][cH]1",
            "fragment2": "[cH]1[cH][n][cH][cH][n]1",
            "percentmcs": 0.6666666666666666,
            "radius": 2,
            "smiles1": "c1ccccc1",
            "smiles2": "c1cnccn1",
            "smirks": "[#6:6]1(:[#6:5](:[#6:4](:[#6:3](:[#6:2](:[#6:1]:1-[H])-[H])-[H])-[H])-[H])-[H]>>[#6:6]1(:[#6:5](:[#7:4]:[#6:3](:[#6:2](:[#7:1]:1)-[H])-[H])-[H])-[H]",
            "valid": true
        },
        {
            "fragment1": "[cH]1[cH][cH][cH][cH][cH]1",
            "fragment2": "[cH]1[cH][n][cH][cH][n]1",
            "percentmcs": 0.6666666666666666,
            "radius": 1,
            "smiles1": "c1ccccc1",
            "smiles2": "c1cnccn1",
            "smirks": "[#6:6]1(:[#6:5](:[#6:4](:[#6:3](:[#6:2](:[#6:1]:1-[H])-[H])-[H])-[H])-[H])-[H]>>[#6:6]1(:[#6:5](:[#7:4]:[#6:3](:[#6:2](:[#7:1]:1)-[H])-[H])-[H])-[H]",
            "valid": true
        }
    ]
]
```
