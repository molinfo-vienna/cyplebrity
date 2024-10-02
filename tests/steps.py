from pytest_bdd import given, parsers, when

from cyplebrity import CyplebrityModel


@given("the CYPlebrity model", target_fixture="predictor")
def cyplebrity_model():
    return CyplebrityModel()


@when(
    parsers.parse("the model generates predictions for the molecule representations"),
    target_fixture="predictions",
)
def predictions(representations, predictor, input_type):
    return predictor.predict(
        representations,
        input_type=input_type,
        output_format="record_list",
    )
