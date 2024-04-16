from pytest_bdd import scenario


@scenario(
    "features/consistent_predictions.feature",
    "Predictions stay consistent with previous versions",
)
def test_consistent_predictions():
    pass
