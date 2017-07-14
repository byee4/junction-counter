from base_package import base_module


def test_add_two_numbers_1():
    a = 1
    b = 0
    assert base_module.add_two_numbers(a, b) == 1


def test_add_two_numbers_2():
    a = 0
    b = 1
    assert base_module.add_two_numbers(a, b) == 1