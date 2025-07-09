# Cyplebrity

## Installation

```sh
# fetch github repository
git clone https://github.com/molinfo-vienna/cyplebrity
cd cyplebrity

# install environment
conda env create -f environment.yml
conda activate cyplebrity

# apply fix for old rdkit versions
curl -sS https://gist.githubusercontent.com/shirte/e1734e51dbc72984b2d918a71b68c25b/raw/ae4afece11980f5d7da9e7668a651abe349c357a/rdkit_installation_fix.sh | bash -s cyplebrity

# install pip dependencies
pip install .

# run cyplebrity
cyplebrity "O=C(N=C(N)N)c1sc2c(c(C#N)ccc2)c1"
```

## Contribute

```sh
# install pip dependencies with test packages
pip install -e .[dev,test]

# run tests
ptw
```

## Contributors

* Wojtek Plonka
* Steffen Hirte
* Axinya Tokareva