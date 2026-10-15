import pydada2
import pydada2._core as core


def test_core_versions():
    assert core.core_version() == "0.1.0"
    assert core.upstream_version() == "1.36.0"


def test_package_version():
    assert pydada2.__version__ == "0.1.0"
