"""Run a WD child off the GUI thread while a modal progress dialog pumps Qt events."""
from PyQt5 import QtCore, QtWidgets


class ProcessWorker(QtCore.QThread):
    def __init__(self, io):
        super().__init__()
        self.io = io
        self.error = None

    def run(self):
        try:
            self.io._run_process()
        except Exception as error:
            self.error = error


class JobDialog(QtWidgets.QDialog):
    def __init__(self, io):
        super().__init__(QtWidgets.QApplication.activeWindow())
        self.io = io
        self.setWindowTitle("WD calculation")
        layout = QtWidgets.QVBoxLayout(self)
        self.label = QtWidgets.QLabel("Calculating...", self)
        layout.addWidget(self.label)
        bar = QtWidgets.QProgressBar(self)
        bar.setRange(0, 0)
        layout.addWidget(bar)
        layout.addWidget(QtWidgets.QLabel("Time limit (minutes):", self))
        limit = QtWidgets.QDoubleSpinBox(self)
        limit.setRange(0.001, 10080)
        limit.setDecimals(3)
        limit.setValue(io.timeout / 60)
        limit.valueChanged.connect(lambda value: setattr(io, "timeout", value * 60))
        layout.addWidget(limit)
        cancel = QtWidgets.QPushButton("Cancel", self)
        cancel.clicked.connect(self.reject)
        layout.addWidget(cancel)

    def reject(self):
        self.io.cancel()
        self.label.setText("Cancelling...")

    def closeEvent(self, event):
        self.reject()
        event.ignore()


def run_with_progress(io):
    dialog = JobDialog(io)
    worker = ProcessWorker(io)
    loop = QtCore.QEventLoop()
    worker.finished.connect(loop.quit)
    worker.finished.connect(dialog.accept)
    worker.start()
    dialog.setModal(True)
    dialog.show()
    if worker.isRunning():
        loop.exec_()
    dialog.hide()
    dialog.deleteLater()
    worker.deleteLater()
    # finished is delivered only once run() returned; no unbounded GUI wait.
    if worker.error is not None:
        raise worker.error
