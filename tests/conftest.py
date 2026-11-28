import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import pytest
from PyQt5 import QtWidgets


@pytest.fixture(scope="session")
def app():
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    yield app


@pytest.fixture(autouse=True)
def messages(monkeypatch):
    from src.helpers import messenger

    monkeypatch.setattr(messenger.Messenger, "show", lambda self: QtWidgets.QMessageBox.Ok)
