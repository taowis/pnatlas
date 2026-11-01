import pytest

from pnatlas import __version__
from pnatlas.cli import main


def test_version(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
    assert __version__ in capsys.readouterr().out


def test_no_args():
    assert main([]) == 0
