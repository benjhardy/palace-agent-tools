from pathlib import Path

import pytest

DATA = Path(__file__).resolve().parent / "data"


@pytest.fixture
def data():
    return DATA
