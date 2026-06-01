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

def test_prepare_model_with_plugin_present(
    host_guest_mmvt_model_input,
    host_guest_mmvt_model
):
    """
    Test that the model was prepared with the plugin present.
    """
    assert "seekrsamplingmethodsplugin" in host_guest_mmvt_model.plugins