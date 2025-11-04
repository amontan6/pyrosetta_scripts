import pytest
from bootcamp_app import identify_secondary_structure_spans

def test_pytest_runs():
    assert 1 + 1 == 2

ss1= "   EEEEE   HHHHHHHH  EEEEE   IGNOR EEEEEE   HHHHHHHHHHH  EEEEE  HHHH   "
expected1= [(4, 8), (12, 19), (22, 26), (36, 41), (45, 55), (58, 62), (65, 68)]

ss2= "HHHHHHH   HHHHHHHHHHHH      HHHHHHHHHHHHEEEEEEEEEEHHHHHHH EEEEHHH "
expected2= [(1, 7), (11, 22), (29, 40), (41, 50), (51, 57), (59, 62), (63, 65)]

ss3= "EEEEEEEEE EEEEEEEE EEEEEEEEE H EEEEE H H H EEEEEEEE"
expected3= [(1,9), (11, 18), (20, 28), (30, 30), (32, 36), (38, 38), (40, 40), (42, 42), (44, 51)]

@pytest.mark.parametrize("test_input,expected", [(ss1, expected1), (ss2, expected2), (ss3, expected3)])
def test_ss_segments(test_input, expected):
    assert identify_secondary_structure_spans(test_input) == expected