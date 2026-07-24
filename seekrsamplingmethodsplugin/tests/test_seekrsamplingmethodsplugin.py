"""
Unit and regression test for the seekrexampleplugin package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import seekrsamplingmethodsplugin


def test_seekrsamplingmethodsplugin_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "seekrsamplingmethodsplugin" in sys.modules
