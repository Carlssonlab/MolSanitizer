import sys
try:
    from PyQt5.QtWidgets import QApplication, QWidget, QLabel, QLineEdit, QPushButton, QVBoxLayout, QHBoxLayout, QCheckBox, QFileDialog, QMessageBox, QGridLayout
except ImportError as e:
    print("PyQt5 is not installed. Please install it using 'pip install PyQt5'.")
    raise e
import argparse
from MolSanitizer.molSanitizer import Sanitycheck, cleanData

class MolSanitizerGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.setWindowTitle("MolSanitizer GUI")
        self.setGeometry(100, 100, 440, 230)

        # Main layout
        main_layout = QVBoxLayout()

        # Grid layout for main options
        grid_layout = QGridLayout()

        # File selection
        self.file_label = QLabel("Input SMILES file:")
        self.file_entry = QLineEdit(self)
        self.file_button = QPushButton("Browse", self)
        self.file_button.clicked.connect(self.browse_file)

        grid_layout.addWidget(self.file_label, 0, 0)
        grid_layout.addWidget(self.file_entry, 0, 1, 1, 2)
        grid_layout.addWidget(self.file_button, 0, 3)

        # Options
        self.salts_check = QCheckBox("Remove Salts", self)
        self.tautomers_check = QCheckBox("Tautomers", self)
        self.pains_check = QCheckBox("Filter PAINS", self)
        self.stereoisomers_check = QCheckBox("Stereoisomers", self)
        self.protonation_check = QCheckBox("Protonation", self)
        self.db2_check = QCheckBox("Generate DB2", self)
        self.db2_check.stateChanged.connect(self.toggle_db2_options)

        grid_layout.addWidget(self.salts_check, 1, 0)
        grid_layout.addWidget(self.tautomers_check, 1, 1)
        grid_layout.addWidget(self.pains_check, 1, 2)
        grid_layout.addWidget(self.stereoisomers_check, 2, 0)
        grid_layout.addWidget(self.protonation_check, 2, 1)
        grid_layout.addWidget(self.db2_check, 2, 2)

        # DB2 Options (Initially Hidden)
        self.rmsd_label = QLabel("RMSD:")
        self.rmsd_entry = QLineEdit(self)
        self.rmsd_entry.setText("0.25")

        self.rs_label = QLabel("Random seed:")
        self.rs_entry = QLineEdit(self)
        self.rs_entry.setText("42")

        self.numconfs_label = QLabel("Max Confs:")
        self.numconfs_entry = QLineEdit(self)
        self.numconfs_entry.setText("2000")

        self.cleanup_label = QLabel("Cleanup:")
        self.cleanup_check = QCheckBox(self)

        # Add to layout but hide initially
        grid_layout.addWidget(self.rmsd_label, 3, 0)
        grid_layout.addWidget(self.rmsd_entry, 3, 1)
        grid_layout.addWidget(self.numconfs_label, 3, 2)
        grid_layout.addWidget(self.numconfs_entry, 3, 3)
        grid_layout.addWidget(self.rs_label, 4, 0)
        grid_layout.addWidget(self.rs_entry, 4, 1)
        grid_layout.addWidget(self.cleanup_label, 4, 2)
        grid_layout.addWidget(self.cleanup_check, 4, 3)

        self.rmsd_label.hide()
        self.rmsd_entry.hide()
        self.rs_label.hide()
        self.rs_entry.hide()
        self.numconfs_label.hide()
        self.numconfs_entry.hide()
        self.cleanup_label.hide()
        self.cleanup_check.hide()

        # Run button at the bottom
        self.run_button = QPushButton("Run MolSanitizer", self)
        self.run_button.clicked.connect(self.run_molSanitizer)
        run_button_layout = QHBoxLayout()
        run_button_layout.addStretch(1)
        run_button_layout.addWidget(self.run_button)
        run_button_layout.addStretch(1)

        # Add layouts to the main layout
        main_layout.addLayout(grid_layout)
        main_layout.addLayout(run_button_layout)

        self.setLayout(main_layout)

    def browse_file(self):
        options = QFileDialog.Options()
        file_path, _ = QFileDialog.getOpenFileName(self, "Select SMILES File", "", "SMILES Files (*.txt *.smi *.ism);;All Files (*)", options=options)
        if file_path:
            self.file_entry.setText(file_path)

    def toggle_db2_options(self):
        if self.db2_check.isChecked():
            self.rmsd_label.show()
            self.rmsd_entry.show()
            self.rs_label.show()
            self.rs_entry.show()
            self.numconfs_label.show()
            self.numconfs_entry.show()
            self.cleanup_label.show()
            self.cleanup_check.show()
        else:
            self.rmsd_label.hide()
            self.rmsd_entry.hide()
            self.rs_label.hide()
            self.rs_entry.hide()
            self.numconfs_label.hide()
            self.numconfs_entry.hide()
            self.cleanup_label.hide()
            self.cleanup_check.hide()

    def run_molSanitizer(self):
        input_file = self.file_entry.text()
        if not input_file:
            QMessageBox.critical(self, "Input Error", "Please select a SMILES file.")
            return

        args = argparse.Namespace(
            input_files=[input_file],
            smiles=None,
            removesalts=self.salts_check.isChecked(),
            tautomers=self.tautomers_check.isChecked(),
            pains=self.pains_check.isChecked(),
            stereoisomers=self.stereoisomers_check.isChecked(),
            protonation=self.protonation_check.isChecked(),
            enamine=False,
            debug=False,
            db2=self.db2_check.isChecked(),
            prefix=None,
            test=False,
            create_custom=False,
            custom=None,
            unwanted=None,
            max_isomers=None,
            randomSeed=int(self.rs_entry.text()) if self.db2_check.isChecked() else None,
            numconfs=int(self.numconfs_entry.text()) if self.db2_check.isChecked() else None,
            rmsd=float(self.rmsd_entry.text()) if self.db2_check.isChecked() else None,
            cleanup=self.cleanup_check.isChecked(),
            lazy=False
        )

        # Run the existing MolSanitizer logic
        try:
            args = Sanitycheck(args)
            cleanData(args)
            QMessageBox.information(self, "Success", "MolSanitizer completed successfully!")
        except Exception as e:
            QMessageBox.critical(self, "Execution Error", str(e))


def main():
    app = QApplication(sys.argv)
    ex = MolSanitizerGUI()
    ex.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
