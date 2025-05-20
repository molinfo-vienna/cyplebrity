@cyplebrity
Feature: Valid predictions

  Scenario Outline: Predictions are valid
    Given a random seed set to <seed>
    And a list of <num_molecules> random molecules, where <num_none> entries are None
    And the representations of the molecules in <input_type> format
    And the CYPlebrity model

    When the model generates predictions for the molecule representations
    And the subset of the result where the input was not None is considered

    Then The result should contain the columns:
            prediction_1
            neighbor_1
            prediction_2
            neighbor_2
            prediction_3
            neighbor_3
            prediction_4
            neighbor_4
            prediction_5
            neighbor_5
    # TODO: filter invalid molecules
    # And the value in column 'prediction_1' should be between 0 and 1
    # And the value in column 'neighbor_1' should be between 0 and 1
    # And the value in column 'prediction_2' should be between 0 and 1
    # And the value in column 'neighbor_2' should be between 0 and 1
    # And the value in column 'prediction_3' should be between 0 and 1
    # And the value in column 'neighbor_3' should be between 0 and 1
    # And the value in column 'prediction_4' should be between 0 and 1
    # And the value in column 'neighbor_4' should be between 0 and 1
    # And the value in column 'prediction_5' should be between 0 and 1
    # And the value in column 'neighbor_5' should be between 0 and 1


  Examples:
  | seed | num_molecules | num_none | input_type |
  | 1    | 10            | 0        | smiles     |
  | 2    | 10            | 1        | smiles     |
  | 3    | 10            | 2        | smiles     |
  | 4    | 10            | 9        | smiles     |
  | 1    | 10            | 0        | smiles     |
  | 2    | 10            | 1        | smiles     |
  | 3    | 10            | 2        | smiles     |
  | 4    | 10            | 9        | smiles     |
