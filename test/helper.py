import pytest
import os
import shlex
import subprocess
import sys
import tempfile


MODULE = "pairef"


def run(line, **kwargs):
    # try:
    #      if os.environ["TRAVIS"] == "true":
    #         # line += (" -t")
    #         # print("Light-testing for TRAVIS CI - not executing REFMAC5")
    #         pass
    # except KeyError:
    #     # print("Normal testing mode.")
    #     pass
    print('\n$ cctbx.python -m pairef', line)
    # print('\n$ python -m', MODULE, line)
    command = ["cctbx.python", '-m', "pairef"] + shlex.split(line)
    # return subprocess.run(command,
    process = subprocess.Popen(command,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True,
                               **kwargs)
    cp = AttrDict()
    cp.stdout, cp.stderr = process.communicate()
    cp.returncode = process.returncode
    return cp


def run_ok(*args, **kwargs):
    cp = run(*args, **kwargs)
    assert cp.returncode == 0
    assert not cp.stderr
    print(cp.stdout)
    return cp


def config(name):
    return os.path.dirname(os.path.realpath(__file__)) + "/fixtures/" + name
    # return pathlib.Path(__file__).parent / 'fixtures' / name


@pytest.fixture
def tmp_environ():
    os.chdir(tempfile.gettempdir())
    return


# https://stackoverflow.com/
# questions/4984647/accessing-dict-keys-like-an-attribute
class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self
