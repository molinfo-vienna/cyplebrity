from pytest_bdd import given, parsers

from cyplebrity import CyplebrityModel


@given("the CYPlebrity model", target_fixture="model")
def model():
    return CyplebrityModel()


@given(parsers.parse("the input type is '{input_type}'"), target_fixture="input_type")
def input_type(input_type):
    return input_type
