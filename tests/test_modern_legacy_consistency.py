import pytest
import os
import shutil
import tempfile
from pathlib import Path
from datpy.datpy import (
    run_legacy_datp,
    run_datp,
)


def test_modern_legacy_consistency():
    thisfile_dir = Path(__file__).parent.resolve()
    inp_path = (thisfile_dir.parent / 'legacy-tests' / 'test_001' / 'input').resolve()
    with tempfile.TemporaryDirectory() as tmpdirname:
        print(tmpdirname)
        shutil.copy(inp_path / 'DAT.INP', tmpdirname) 
        shutil.copy(inp_path / 'GMDATA.CRD', tmpdirname) 
        os.chdir(tmpdirname)
        run_legacy_datp(dbfile_out='input.json', do_reduce=False)
        run_legacy_datp(dbfile_out='legacy_result.json', do_reduce=True)
        run_datp('input.json', 'modern_result.json')
        with open('modern_result.json', 'r') as f:
            modern_result_text = f.read()
        with open('legacy_result.json', 'r') as f:
            legacy_result_text = f.read()
    assert modern_result_text == legacy_result_text 
