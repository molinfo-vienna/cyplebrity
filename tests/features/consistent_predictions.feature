@cyplebrity
Feature: Consistent predictions

  Scenario Outline: Predictions stay consistent with previous versions
    Given an input molecule specified by '<input_smiles>'
    And the CYPlebrity model
    
    When the model generates predictions for the molecule representations
    And the subset of the result where the input was not None is considered
    
    Then the value in column 'name' should be equal to '<name>'
    And the value in column 'input_text' should be equal to '<input_smiles>'
    And the value in column 'prediction_1' should be equal to <prediction_1>
    And the value in column 'neighbor_1' should be equal to <neighbor_1>
    And the value in column 'prediction_2' should be equal to <prediction_2>
    And the value in column 'neighbor_2' should be equal to <neighbor_2>
    And the value in column 'prediction_3' should be equal to <prediction_3>
    And the value in column 'neighbor_3' should be equal to <neighbor_3>
    And the value in column 'prediction_4' should be equal to <prediction_4>
    And the value in column 'neighbor_4' should be equal to <neighbor_4>
    And the value in column 'prediction_5' should be equal to <prediction_5>
    And the value in column 'neighbor_5' should be equal to <neighbor_5>

    Examples:
    | name                     | input_smiles                                                                                                  | preprocessed_smiles                                                                             | prediction_1 | neighbor_1 | prediction_2 | neighbor_2 | prediction_3 | neighbor_3 | prediction_4 | neighbor_4 | prediction_5 | neighbor_5 | errors |
    | Aciclovir                | C1=NC2=C(N1COCCO)N=C(NC2=O)N Aciclovir                                                                        | N=c1[nH]c(=O)c2ncn(COCCO)c2[nH]1                                                                | 0.21         | 0.37       | 0.07         | 0          | 0.05         | 0          | 0.06         | 0          | 0.06         | 0          |        |
    | Amiodarone               | CCN(CC)CCOc1c(I)cc(cc1I)C(=O)c2c3ccccc3oc2CCCC Amiodarone                                                     | CCCCc1oc2ccccc2c1C(=O)c1cc(I)c(OCCN(CC)CC)c(I)c1                                                | 0.71         | 0.28       | 0.87         | 0.28       | 0.48         | 0          | 0.62         | 0.28       | 0.56         | 0.28       |        |
    | Arsphenamine (Salvarsan) | C1=CC(=C(C=C1[As]=[As]C2=CC(=C(C=C2)O)N)N)O.Cl.Cl Arsphenamine (Salvarsan)                                    | None                                                                                            | None         | None       | None         | None       | None         | None       | None         | None       | None         | None       | !1     |
    | Cyclophosphamide         | C1CNP(=O)(OC1)N(CCCl)CCCl Cyclophosphamide                                                                    | O=P1(N(CCCl)CCCl)NCCCO1                                                                         | 0.09         | 0          | 0.09         | 0          | 0.09         | 0          | 0.14         | 0          | 0.07         | 0          |        |
    | Doxorubicin              | C[C@H]1[C@H]([C@H](C[C@@H](O1)O[C@H]2C[C@@](Cc3c2c(c4c(c3O)C(=O)c5cccc(c5C4=O)OC)O)(C(=O)CO)O)N)O Doxorubicin | COc1cccc2c(O)c3c(O)c4c(c(O)c3c(O)c12)=C(O[C@H]1C[C@H](N)[C@H](O)[C@H](C)O1)C[C@](O)(C(=O)CO)C=4 | 0.19         | 0          | 0.05         | 0          | 0.07         | 0          | 0.18         | 0.23       | 0.2          | 0.18       |        |
    | Hydrochlorothiazide      | O=S(=O)(N)c1c(Cl)cc2c(c1)S(=O)(=O)NCN2 Hydrochlorothiazide                                                    | NS(=O)(=O)c1cc2c(cc1Cl)NCNS2(=O)=O                                                              | 0.1          | 0          | 0.13         | 0          | 0.1          | 0          | 0.12         | 0          | 0.07         | 0          |        |
    | Levofloxacin             | C[C@H]1COc2c3n1cc(c(=O)c3cc(c2N4CCN(CC4)C)F)C(=O)O Levofloxacin                                               | C[C@H]1COc2c(N3CCN(C)CC3)c(F)cc3c(=O)c(C(=O)O)cn1c23                                            | 0.63         | 0.24       | 0.22         | 0.54       | 0.08         | 0          | 0.09         | 0          | 0.21         | 0.54       |        |
    | Metamizole (Sulpyrine)   | CC1=C(C(=O)N(N1C)C2=CC=CC=C2)N(C)CS(=O)(=O)[O-].O.[Na+] Metamizole (Sulpyrine)                                | Cc1c(N(C)CS(=O)(=O)[O-])c(=O)n(-c2ccccc2)n1C                                                    | 0.17         | 0          | 0.17         | 0          | 0.15         | 0          | 0.1          | 0          | 0.14         | 0          |        |
    | Nifedipine               | CC1=C(C(C(=C(N1C=O)C)C(=O)OC)C2=CC=CC=C2[N+](=O)[O-])C(=O)OC Nifedipine                                       | COC(=O)C1=C(C)N(C=O)C(C)=C(C(=O)OC)C1c1ccccc1[N+](=O)[O-]                                       | 0.53         | 0.6        | 0.36         | 0.64       | 0.48         | 0.64       | 0.19         | 0.64       | 0.31         | 0.64       |        |
    | Phenprocoumon            | OC=1c3ccccc3OC(=O)C=1C(CC)c2ccccc2 Phenprocoumon                                                              | CCC(c1ccccc1)c1c(O)c2ccccc2oc1=O                                                                | 0.39         | 0.29       | 0.94         | 0          | 0.31         | 0.37       | 0.15         | 0          | 0.19         | 0.29       |        |
    | Rivaroxaban              | O=C1COCCN1c2ccc(cc2)N3C[C@@H](OC3=O)CNC(=O)c4ccc(s4)Cl Rivaroxaban                                            | O=C(NC[C@H]1CN(c2ccc(N3CCOCC3=O)cc2)C(=O)O1)c1ccc(Cl)s1                                         | 0.41         | 0.38       | 0.29         | 0.38       | 0.3          | 0.38       | 0.45         | 0.38       | 0.59         | 0.22       |        |