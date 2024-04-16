from pytest_bdd import scenario


@scenario(
    "features/valid_predictions.feature",
    "Predictions are valid",
)
def test_valid_predictions():
    pass
