import pytest
from bootcamp_app import identify_secondary_structure_spans, fold_tree_from_dssp_string

def test_pytest_runs():
    assert 1 + 1 == 2

ss1= "   EEEEE   HHHHHHHH  EEEEE   IGNOR EEEEEE   HHHHHHHHHHH  EEEEE  HHHH   "
expected1= [(4, 8), (12, 19), (22, 26), (36, 41), (45, 55), (58, 62), (65, 68)]

ss2= "HHHHHHH   HHHHHHHHHHHH      HHHHHHHHHHHHEEEEEEEEEEHHHHHHH EEEEHHH "
expected2= [(1, 7), (11, 22), (29, 40), (41, 50), (51, 57), (59, 62), (63, 65)]

ss3= "EEEEEEEEE EEEEEEEE EEEEEEEEE H EEEEE H H H EEEEEEEE"
expected3= [(1,9), (11, 18), (20, 28), (30, 30), (32, 36), (38, 38), (40, 40), (42, 42), (44, 51)]

ft_input = "   EEEEEEE    EEEEEEE         EEEEEEEEE    EEEEEEEEEE   HHHHHH         EEEEEEEEE         EEEEE     "
ft_expected = [(7, 1, -1), (7, 10, -1), (7, 12, 1), (12, 11, -1), (12, 14, -1), (7, 18, 2), (18, 15, -1), (18, 21, -1), (7, 26, 3), (26, 22, -1), (26, 30, -1), (7, 35, 4), (35, 31, -1), (35, 39, -1), (7, 41, 5), (41, 40, -1), (41, 43, -1), (7, 48, 6), (48, 44, -1), (48, 53, -1), (7, 55, 7), (55, 54, -1), (55, 56, -1), (7, 59, 8), (59, 57, -1), (59, 62, -1), (7, 67, 9), (67, 63, -1), (67, 71, -1), (7, 76, 10), (76, 72, -1), (76, 80, -1), (7, 85, 11), (85, 81, -1), (85, 89, -1), (7, 92, 12), (92, 90, -1), (92, 99, -1)]

@pytest.mark.parametrize("test_input,expected", [(ss1, expected1), (ss2, expected2), (ss3, expected3)])
def test_ss_segments(test_input, expected):
    assert identify_secondary_structure_spans(test_input) == expected

@pytest.mark.parametrize("test_input,expected", [(ft_input, ft_expected)])
def test_fold_tree_from_string(test_input, test_expected):
    assert fold_tree_from_dssp_string(test_input) == test_expected