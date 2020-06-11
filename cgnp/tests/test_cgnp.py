"""
Unit and regression test for the cgnp package.
"""

# Import package, test suite, and other packages as needed
import cgnp
import pytest
import sys
import mbuild as mb

def test_cgnp_imported():
    """ Sample test, will always pass so long as import statement worked """
    assert "cgnp" in sys.modules

def test_import():
    """ Test that mBuild recipe import works """
    assert "cgnp" in vars(mb.recipes).keys()
