from pytest_bdd import given

from cyplebrity import CyplebrityModel


@given("the CYPlebrity model", target_fixture="model")
def cyplebrity_model():
    return CyplebrityModel()
