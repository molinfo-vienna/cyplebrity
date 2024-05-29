from nerdd_module import auto_cli

from .cyplebrity_model import CyplebrityModel


@auto_cli
def main():
    return CyplebrityModel()
