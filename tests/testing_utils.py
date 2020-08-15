import os
from pkg_resources import Requirement, resource_filename, ResolutionError


def get_test_data(filename):
    filepath = None
    try:
        filepath = resource_filename(Requirement.parse("pyfastx"), "tests/data/" + filename)
    except ResolutionError:
        pass
    if not filepath or not os.path.isfile(filepath):
        filepath = os.path.join(os.path.dirname(__file__), 'data', filename)
    return filepath
